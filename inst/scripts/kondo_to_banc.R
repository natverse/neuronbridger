#!/usr/bin/env Rscript
#
# kondo_to_banc.R --- HISTORICAL RECORD (not a user-facing tool)
#
# Bridges the Kondo et al. 2020 endogenous neurotransmitter-receptor
# light-microscopy stacks from IS2 -> FCWB -> JRC2018F -> BANC voxel
# space and writes each as a Neuroglancer "precomputed" image layer
# named `{gene}_{sample}.ng/`. This file is the source of record for
# the LM volumes hosted under
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
#                                       kondo_et_al_2020/
# and indexed by the registry at
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
#                                       registry.json
# Browse the registry from R via `bancr:::banc_lm_volumes()`.
#
# It is intentionally NOT exposed via roxygen / package docs --- end
# users do not need to re-run it. Re-run is only relevant if the
# upstream data, the registration chain, or the lm_to_banc_layer()
# orchestrator changes.
#
# Usage (lee-lab maintainers, on a machine that has the upstream
# NRRDs cached locally + transformix on PATH + the Saalfeld
# template-building JAR):
#
#   Rscript inst/scripts/kondo_to_banc.R \
#     [test|glutamate|all] [--upload] [--register] [--data-root <path>]
#
#     test            --  3 volumes  --- runtime probe
#     glutamate       -- 22 volumes  --- 11 glutamate-receptor genes x 2 samples
#     all             -- 142 volumes --- everything in the Kondo inventory
#     --upload        --  on completion, gsutil cp -Z .ng dirs to GCS
#                         (gs://lee-lab.../light_level/kondo_et_al_2020/).
#     --register      --  on completion, merge registry_entries.json into
#                         the master registry.json on GCS, then re-mint
#                         bancr::banc_lm_links via the bancr data-raw helper.
#                         Implies --upload.
#     --data-root P   --  override the local data root. Defaults to
#                         ~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat
#                         (or $NEURONBRIDGER_DATA_ROOT if set).
#     --one <name>    --  worker mode: process one sample "<gene>_<sample>"
#                         (e.g. "GluRIIA_no1"), write per-sample RDS
#                         stubs, exit. Used internally by outer mode so
#                         each sample runs in a fresh R heap.
#
# Local layout expected under <data-root>:
#   kondo_et_al_2020/nrrd/    G-Node IS2-aligned NRRDs (sources)
#   templates/                JRC2018F.nrrd etc.
#   kondo_to_banc/            outputs (this script writes <DATASET>/ here)
#
# The script is resumable: any `<gene>_<sample>.ng/` directory that
# already exists is silently skipped.
#
# == Data source ==========================================================
#
# Kondo, S., Takahashi, T., Yamagata, N., Imanishi, Y., Katow, H.,
# Hiramatsu, S., Lynn, K., Abe, A., Kumaraswamy, A., Tanimoto, H.
# (2020). Neurochemical organization of the Drosophila brain visualized
# by endogenously tagged neurotransmitter receptors. Cell Rep
# 30(1):284-297.e5. PMID 31914394.
#
# IS2-aligned receptor NRRDs (one stack per gene x specimen x channel)
# are hosted on G-Node:
#   https://doi.gin.g-node.org/10.12751/g-node.10246f/
#
# Filename convention (in the G-Node distribution):
#   IS2_<gene>_no<sample>_<channel>_warp_m0g40c4e1e-1x16r3.nrrd
#     <sample>  is `no1` or `no2`           --- the two prepared specimens
#     <channel> is `01` (NC82 reference) or `02` (GFP / receptor)
#
# We always pull channel `02` (the GFP-tagged receptor signal). Channel
# `01` is the NC82 anatomical reference and is not mirrored to BANC.
#
# == Pipeline ============================================================
#
# Per-volume image-mode warp chain (see ?neuronbridger::lm_pipeline):
#   Stage A  IS2  -> FCWB        --- CMTK reformatx (nat.flybrains)
#   Stage B  FCWB -> JRC2018F    --- saalfeldlab template-building
#                                    RenderTransformed (H5 displacement
#                                    field; nat.jrcbrains)
#   Stage C  JRC2018F -> BANC    --- Elastix transformix (bancr's
#                                    bundled brain_240721 chain)
#   Stage D  BANC NRRD -> precomputed --- nrrd_to_precomputed()
#
# Per-volume runtime on Apple-silicon Mac with 8 threads is ~13 min
# (Stage A 7 min, B 3 min, C 2 min, D 30 s); CMTK Stage A is the
# single-threaded bottleneck.
#
# Output written to: OUT_ROOT/<dataset>/<gene>_<sample>.ng/  +
# per-batch timings.csv and registry_entries.json. Upload to GCS and
# registry-merge are out-of-band; the script prints the gsutil/jq
# commands at the end.
#
# == Outputs (cumulative across all runs) ================================
#
# 22 glutamate-receptor layers (May 2026 batch); the rest of the
# inventory is being processed in follow-up runs. Each layer is on the
# BANC voxel grid (2400 x 924 x 789 @ 400 nm), uint8, gzip-compressed
# 64^3 chunks. Median fg-voxel inside-mesh fraction (vs.
# bancr::banc_brain_neuropil.surf) ~75-97% across the glutamate set.

suppressMessages({
  library(neuronbridger)
  library(nat.jrcbrains)
  library(jsonlite)
})

# --- knobs ------------------------------------------------------------

# All local sources + outputs live under a single DATA_ROOT.
DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"
DATASET    <- "kondo_et_al_2020"
GS_BASE    <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level"
BANCR_REPO <- "~/Projects/flyconnectome/bancr"

# --- args -----------------------------------------------------------
.args  <- commandArgs(trailingOnly = TRUE)
dr_idx <- which(.args == "--data-root")
data_root <- if (length(dr_idx) && dr_idx < length(.args)) .args[dr_idx + 1L] else Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT)
DATA_ROOT <- path.expand(data_root)
if (length(dr_idx)) .args <- .args[-c(dr_idx, dr_idx + 1L)]

force    <- "--force" %in% .args
upload   <- "--upload"   %in% .args || "--register" %in% .args
register <- "--register" %in% .args
one_idx  <- which(.args == "--one")
one      <- if (length(one_idx) && one_idx < length(.args)) .args[one_idx + 1L] else NA_character_
if (length(one_idx)) .args <- .args[-c(one_idx, one_idx + 1L)]
mode <- {
  m <- .args[!.args %in% c("--force", "--upload", "--register")]
  if (!length(m)) "test" else m[1]
}

DROP_DIR <- file.path(DATA_ROOT, "kondo_et_al_2020", "nrrd")
OUT_ROOT <- file.path(DATA_ROOT, "kondo_to_banc")
cat("DATA_ROOT: ", DATA_ROOT, "\n", sep = "")

# --- gene cohorts (Kondo 2020 inventory, by neurotransmitter family) --

GENES_GLUTAMATE <- c(
  # ionotropic AMPA-like (NMJ-pre/post receptors)
  "GluRIIA","GluRIIB","GluRIIC","GluRIID","GluRIIE","GluRIB",
  # NMDA-type
  "Nmdar1","Nmdar2",
  # ionotropic chloride-permeable (glutamate-gated chloride)
  "GluClalpha",
  # metabotropic
  "mGluR",
  # kainate
  "DKaiR1D"
)

GENES_ACETYLCHOLINE <- c(
  # nicotinic alpha
  "nAChRa1","nAChRa2","nAChRa3","nAChRa6","nAChRa7",
  # nicotinic beta
  "nAChRb1","nAChRb2","nAChRb3",
  # muscarinic
  "mAChR-A","mAChR-B"
)

GENES_GABA_HISTAMINE_GLYCINE <- c(
  # GABA-B metabotropic
  "GABA-B-R1","GABA-B-R2",
  # GABA-A / glycine-like ionotropic (GluCl-family)
  "Grd","Lcch3",
  # histamine-gated chloride
  "ort","hec",
  # other
  "pdfr"
)

GENES_BIOGENIC_AMINES <- c(
  # serotonin
  "5-HT1B","5-HT2A","5-HT2B","5-HT7",
  # dopamine
  "DopR1","DopR2","DopEcR","D2R",
  # octopamine / tyramine
  "Oamb","Octb3R","TyrRII"
)

GENES_NEUROPEPTIDES <- c(
  "AkhR",
  "AstA-R1","AstA-R2","AstC-R1","AstC-R2",
  "CapaR","CCAP-R","CCHa1-R","CCHa2-R",
  "CCKLR-17D1","CCKLR-17D3",
  "CrzR",
  "Dh31-R","Dh44-R2",
  "NPFR","PK1-R","PK2-R2","RYa-R",
  "SIFaR","SPR","TkR86C"
)

GENES_OTHER <- c(
  # leucine-rich-repeat GPCRs
  "lgr1","Lgr3","lgr4",
  # unclassified / orphan
  "rk","mtt",
  # CG-numbered / uncharacterised
  "CG7589","CG12344","CG13229","CG13575","CG13995","CG32547","CG43795"
)

GENES_ALL <- c(GENES_GLUTAMATE, GENES_ACETYLCHOLINE,
               GENES_GABA_HISTAMINE_GLYCINE, GENES_BIOGENIC_AMINES,
               GENES_NEUROPEPTIDES, GENES_OTHER)

# Genes for which the upstream Kondo distribution has only sample no1
# (no2 absent on G-Node). They contribute one volume to the cohort.
GENES_SINGLE_SAMPLE <- c("D2R","PK1-R")

# Kondo file-naming sometimes uses synonyms; map the user-facing
# registry name (LHS) to the actual G-Node filename stem (RHS).
GENE_FILE_SYNONYM <- list(
  DKaiR1D    = "KaiR1D",
  GluClalpha = "GluCla"
)

# As of 2026-05-05 the local mirror has a number of corrupt
# (gzip-truncated, crc-mismatched) NRRDs. They cannot be warped
# through CMTK / RenderTransformed without a re-download from G-Node.
# Listed here as `(gene, sample)` pairs; the script skips them with a
# CORRUPT message so the rest of the cohort can still be processed.
CORRUPT_GENE_SAMPLE <- list(
  c("nAChRa1", "no1"), c("nAChRa1", "no2"),
  c("nAChRa2", "no1"), c("nAChRa2", "no2"),
  c("nAChRa3", "no1"), c("nAChRa3", "no2"),
  c("nAChRa6", "no1"), c("nAChRa6", "no2"),
  c("nAChRa7", "no1"), c("nAChRa7", "no2"),
  c("nAChRb1", "no1"), c("nAChRb1", "no2"),
  c("nAChRb2", "no1"), c("nAChRb2", "no2"),
  c("nAChRb3", "no1"), c("nAChRb3", "no2"),
  c("Dh31-R",  "no1"), c("Dh31-R",  "no2"),
  c("Octb3R",  "no1"), c("Octb3R",  "no2"),
  c("Lgr3",    "no1"),
  c("AstA-R2", "no2")
)
.is_corrupt <- function(gene, sample) {
  any(vapply(CORRUPT_GENE_SAMPLE, function(p)
             identical(p[[1]], gene) && identical(p[[2]], sample),
             logical(1)))
}

# Build per-(gene, sample) target list given a vector of gene names.
.expand_targets <- function(genes) {
  out <- list()
  for (g in genes) {
    samples <- if (g %in% GENES_SINGLE_SAMPLE) "no1" else c("no1","no2")
    for (s in samples) out[[length(out) + 1L]] <- c(g, s)
  }
  out
}

GENE_SET_TEST       <- list(c("GluRIIA","no1"), c("Nmdar1","no1"), c("mGluR","no1"))
GENE_SET_GLUTAMATE  <- .expand_targets(GENES_GLUTAMATE)
GENE_SET_ALL        <- .expand_targets(GENES_ALL)

target_set <- switch(mode[1],
                     test      = GENE_SET_TEST,
                     glutamate = GENE_SET_GLUTAMATE,
                     all       = GENE_SET_ALL,
                     stop("mode must be 'test', 'glutamate', or 'all'; got: ",
                          mode[1]))

# --- locate sources ---------------------------------------------------

find_source <- function(gene, sample) {
  candidate_names <- unique(c(gene, GENE_FILE_SYNONYM[[gene]]))
  for (nm in candidate_names) {
    pat <- sprintf("(^|_)%s_%s_02_warp_.*\\.nrrd$", nm, sample)
    cand <- list.files(DROP_DIR, pattern = pat, full.names = TRUE)
    if (length(cand)) return(cand[1])
  }
  NA_character_
}

if (!dir.exists(DROP_DIR))
  stop("Kondo NRRD dir not found: ", DROP_DIR)

out_root <- normalizePath(path.expand(OUT_ROOT), mustWork = FALSE)
out_dir  <- file.path(out_root, DATASET)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stubs_dir <- file.path(out_dir, "_per_sample")
dir.create(stubs_dir, showWarnings = FALSE, recursive = TRUE)

# --- worker mode: one sample, fresh R heap, drop RDS stubs, exit ---
if (!is.na(one)) {
  # `one` is "<gene>_<sample>" (e.g. "GluRIIA_no1"). Split on the LAST
  # underscore so gene names containing "_" survive; Kondo genes don't
  # normally, but a "no1"/"no2" tail is deterministic.
  parts <- regmatches(one, regexec("^(.*)_(no\\d+)$", one))[[1]]
  if (length(parts) < 3L) stop("--one must be '<gene>_<sample>', got: ", one)
  gene <- parts[[2]]; sample <- parts[[3]]
  name <- paste(gene, sample, sep = "_")
  if (.is_corrupt(gene, sample)) stop("source marked corrupt: ", name)
  src <- find_source(gene, sample)
  if (is.na(src)) stop("no source NRRD for ", name)
  nat.jrcbrains::register_saalfeldlab_registrations()
  message("--- ", name, " (worker, PID=", Sys.getpid(), ") ---")
  res <- lm_to_banc_layer(
    input        = src,
    gene         = gene,
    sample       = sample,
    channel      = "02",
    dataset      = DATASET,
    output_dir   = out_dir,
    source_path  = src
  )
  saveRDS(data.frame(name = name,
                     stageA = res$timings[["stageA_is2_to_fcwb"]],
                     stageB = res$timings[["stageB_fcwb_to_jrcf"]],
                     stageC = res$timings[["stageC_jrcf_to_banc"]],
                     stageD = res$timings[["stageD_precomputed"]],
                     total  = res$timings[["total"]]),
          file.path(stubs_dir, paste0(name, ".timings.rds")))
  saveRDS(res$registry_entry,
          file.path(stubs_dir, paste0(name, ".registry.rds")))
  quit(save = "no", status = 0L)
}

# --- run --------------------------------------------------------------

processed_names <- character(0); failed_names <- character(0)
self <- sub("^--file=", "",
            commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1L])
self <- normalizePath(self, mustWork = TRUE)

for (gs in target_set) {
  gene   <- gs[[1]]; sample <- gs[[2]]
  name   <- paste(gene, sample, sep = "_")
  if (.is_corrupt(gene, sample)) {
    message(sprintf("SKIP  %s --- source file is corrupt (re-download needed)", name))
    next
  }
  src    <- find_source(gene, sample)
  if (is.na(src)) {
    message(sprintf("SKIP  %s --- no source file found", name))
    next
  }
  pc_dir <- file.path(out_dir, paste0(name, ".ng"))
  if (dir.exists(pc_dir) && length(list.files(pc_dir)) && !force) {
    message(sprintf("SKIP  %s (precomputed dir exists; pass --force to recompute)", name))
    if (file.exists(file.path(stubs_dir, paste0(name, ".registry.rds"))))
      processed_names <- c(processed_names, name)
    next
  }
  worker_args <- c(shQuote(self), "--one", shQuote(name),
                   "--data-root", shQuote(DATA_ROOT))
  if (force) worker_args <- c(worker_args, "--force")
  rc <- system2("Rscript", worker_args)
  if (rc != 0L) {
    message("  FAILED: ", name, " (worker exit ", rc, ")")
    failed_names <- c(failed_names, name); next
  }
  processed_names <- c(processed_names, name)
}

if (length(failed_names))
  cat(sprintf("\n%d sample(s) FAILED: %s\n",
              length(failed_names), paste(failed_names, collapse = ", ")))

# Aggregate per-sample RDS stubs into batch timings + registry
timings <- lapply(processed_names, function(nm) {
  f <- file.path(stubs_dir, paste0(nm, ".timings.rds"))
  if (file.exists(f)) readRDS(f) else NULL
})
timings <- Filter(Negate(is.null), timings)
results_reg <- lapply(processed_names, function(nm) {
  f <- file.path(stubs_dir, paste0(nm, ".registry.rds"))
  if (file.exists(f)) readRDS(f) else NULL
})
results_reg <- Filter(Negate(is.null), results_reg)

if (!length(results_reg)) {
  message("\nNo new volumes processed (all already on disk).")
  quit(save = "no", status = 0L)
}
if (!length(timings)) { cat("\nNo timings — RDS stubs missing?\n"); quit(save = "no") }

# --- emit per-batch timings + registry-entries stub -------------------

times_df <- do.call(rbind, timings)
write.csv(times_df, file.path(out_dir, "timings.csv"), row.names = FALSE)

cat("\n=== timings (s) ===\n")
print(times_df, row.names = FALSE)

reg_entries <- do.call(rbind, results_reg)
reg_json    <- list(
  schema_version = 1L,
  volumes        = lapply(seq_len(nrow(reg_entries)), function(i) {
    row <- as.list(reg_entries[i, ])
    # voxdims_nm is a list-column; flatten to a plain integer vector
    # so the JSON emits [400,400,400] not [[400,400,400]].
    row$voxdims_nm <- as.integer(unlist(row$voxdims_nm))
    row
  })
)
local_reg <- file.path(out_dir, "registry_entries.json")
writeLines(jsonlite::toJSON(reg_json, auto_unbox = TRUE, pretty = TRUE),
           local_reg)

# --- optional upload + registry merge ------------------------------
if (upload) {
  cat("\n=== uploading .ng dirs to GCS ===\n")
  cmd <- sprintf("gsutil -m cp -Z -r %s/*.ng %s/%s/",
                 shQuote(out_dir), GS_BASE, DATASET)
  cat("$ ", cmd, "\n", sep = "")
  if (system(cmd) != 0L) stop("gsutil cp failed")
}
if (register) {
  cat("\n=== merging into master registry.json ===\n")
  master_url <- sprintf("%s/registry.json", GS_BASE)
  tmp_master <- tempfile(fileext = ".json")
  tmp_merged <- tempfile(fileext = ".json")
  if (system2("gsutil", c("cat", master_url),
              stdout = tmp_master) != 0L) stop("gsutil cat failed")
  cmd <- sprintf("jq '.volumes += $batch.volumes' --argfile batch %s %s > %s",
                 shQuote(local_reg), shQuote(tmp_master), shQuote(tmp_merged))
  if (system(cmd) != 0L) stop("jq merge failed")
  cmd <- sprintf("gsutil cp %s %s", shQuote(tmp_merged), master_url)
  if (system(cmd) != 0L) stop("gsutil push failed")
  cat("master registry updated\n")

  cat("\n=== re-minting bancr::banc_lm_links ===\n")
  bancr_dir <- normalizePath(path.expand(BANCR_REPO), mustWork = FALSE)
  if (dir.exists(bancr_dir)) {
    cmd <- sprintf("cd %s && Rscript data-raw/make_banc_lm_links.R",
                   shQuote(bancr_dir))
    system(cmd)
  } else {
    cat("(skipping banc_lm_links — bancr repo not at ", bancr_dir, ")\n", sep = "")
  }
}
