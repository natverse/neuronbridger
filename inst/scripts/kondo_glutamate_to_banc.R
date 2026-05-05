#!/usr/bin/env Rscript
#
# Batch driver: bridge a curated Kondo et al. 2020 glutamate-receptor
# subset from IS2 -> JRC2018F -> BANC voxel space and write each as a
# Neuroglancer precomputed image layer named `{gene}_{sample}.ng/`.
#
# Usage:
#   Rscript inst/scripts/kondo_glutamate_to_banc.R [test|all]
#     test  -- 3 volumes (GluRIIA_01, Nmdar1_01, mGluR_01) — runtime probe
#     all   -- 22 volumes (11 genes x 2 samples)
#
# Requires: nat.flybrains, nat.jrcbrains, nat.h5reg, nat.templatebrains,
#           bancr (for inst/extdata/brain_240721/), reticulate +
#           cloud-volume, mmand, transformix on PATH.
#
# Output goes to OUT_ROOT/<dataset>/<gene>_<sample>.ng/ + a per-batch
# timings.csv and registry_entries.json.

suppressMessages({
  library(neuronbridger)
  library(nat.jrcbrains)
  library(jsonlite)
})

# --- knobs ------------------------------------------------------------

DROP_DIR  <- "/Users/asbates/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/kondo_et_al_2020/nrrd"
OUT_ROOT  <- "~/kondo_to_banc"
DATASET   <- "kondo_et_al_2020"

# Kondo file naming: IS2_<gene>_no2_<sample>_warp_m0g40c4e1e-1x16r3.nrrd
# (some files lack the IS2_ prefix — the script handles both).
mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "test"

GENE_SET_TEST <- list(
  c("GluRIIA", "01"),
  c("Nmdar1",  "01"),
  c("mGluR",   "01")
)

GENE_SET_ALL <- (function() {
  genes <- c("GluRIIA","GluRIIB","GluRIIC","GluRIID","GluRIIE",
             "GluRIB","Nmdar1","Nmdar2","GluClalpha","mGluR","DKaiR1D")
  out <- list()
  for (g in genes) for (s in c("01","02")) out[[length(out)+1L]] <- c(g, s)
  out
})()

# Kondo file-naming sometimes uses synonyms / alternate spellings; map
# the user-facing registry name to whatever the Dropbox filename
# actually says. Registry entry uses the LHS; source filename keeps
# the RHS for provenance.
GENE_FILE_SYNONYM <- list(
  DKaiR1D    = "KaiR1D",
  GluClalpha = "GluCla"
)

target_set <- switch(mode[1],
                     test = GENE_SET_TEST,
                     all  = GENE_SET_ALL,
                     stop("mode must be 'test' or 'all', got: ", mode[1]))

# --- locate sources ---------------------------------------------------

find_source <- function(gene, sample) {
  candidate_names <- unique(c(gene, GENE_FILE_SYNONYM[[gene]]))
  for (nm in candidate_names) {
    pat <- sprintf("(^|_)%s_no2_%s_warp_.*\\.nrrd$", nm, sample)
    cand <- list.files(DROP_DIR, pattern = pat, full.names = TRUE)
    if (length(cand)) return(cand[1])
  }
  NA_character_
}

if (!dir.exists(DROP_DIR))
  stop("Kondo NRRD dir not found: ", DROP_DIR)

nat.jrcbrains::register_saalfeldlab_registrations()
out_root <- normalizePath(path.expand(OUT_ROOT), mustWork = FALSE)
out_dir  <- file.path(out_root, DATASET)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- run --------------------------------------------------------------

results  <- list()
timings  <- list()
for (gs in target_set) {
  gene   <- gs[[1]]; sample <- gs[[2]]
  name   <- paste(gene, sample, sep = "_")
  src    <- find_source(gene, sample)
  if (is.na(src)) {
    message(sprintf("SKIP  %s — no source file found", name))
    next
  }
  pc_dir <- file.path(out_dir, paste0(name, ".ng"))
  if (dir.exists(pc_dir) && length(list.files(pc_dir))) {
    message(sprintf("SKIP  %s — precomputed dir already exists", name))
    next
  }
  message(sprintf("\n--- %s ---", name))
  message("  source: ", basename(src))
  res <- try(lm_to_banc_layer(
    input        = src,
    gene         = gene,
    sample       = sample,
    channel      = "no2",
    dataset      = DATASET,
    output_dir   = out_dir,
    source_path  = src
  ), silent = TRUE)
  if (inherits(res, "try-error")) {
    message("  FAILED:\n  ", attr(res, "condition")$message)
    next
  }
  results[[name]] <- res
  timings[[length(timings) + 1L]] <- data.frame(
    name   = name,
    stageA = res$timings[["stageA_is2_to_fcwb"]],
    stageB = res$timings[["stageB_fcwb_to_jrcf"]],
    stageC = res$timings[["stageC_jrcf_to_banc"]],
    stageD = res$timings[["stageD_precomputed"]],
    total  = res$timings[["total"]]
  )
}

if (!length(results)) {
  message("\nNo volumes processed.")
  quit(save = "no", status = 1L)
}

# --- emit timings + registry stub ------------------------------------

times_df <- do.call(rbind, timings)
write.csv(times_df, file.path(out_dir, "timings.csv"), row.names = FALSE)

cat("\n=== timings (s) ===\n")
print(times_df, row.names = FALSE)

reg_entries <- do.call(rbind, lapply(results, `[[`, "registry_entry"))
reg_json    <- list(
  schema_version = 1L,
  volumes        = lapply(seq_len(nrow(reg_entries)), function(i) {
    row <- as.list(reg_entries[i, ])
    # voxdims_nm is a list-column; unwrap to a plain integer vector so
    # the JSON serialiser emits [400,400,400] not [[400,400,400]].
    row$voxdims_nm <- as.integer(unlist(row$voxdims_nm))
    row
  })
)
writeLines(jsonlite::toJSON(reg_json, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "registry_entries.json"))

# --- next-step instructions -----------------------------------------

cat("\nNext steps (lee-lab maintainers only — bucket is public-read but ",
    "private-write):\n", sep = "")
cat(sprintf("  gsutil -m cp -r %s/*.ng gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/%s/\n",
            out_dir, DATASET))
cat("\nThen merge the per-batch registry_entries.json into the master ",
    "registry:\n", sep = "")
cat("  gsutil cat gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json | \\\n")
cat(sprintf("    jq '.volumes += $batch.volumes' --argfile batch %s/registry_entries.json | \\\n",
            out_dir))
cat("    gsutil cp - gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json\n")
