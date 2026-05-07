#!/usr/bin/env Rscript
#
# deng_to_banc.R --- HISTORICAL RECORD (not a user-facing tool)
#
# Bridges Deng et al. 2019 light-microscopy receptor / peptide
# knock-in stacks (raw confocal LSM, native acquisition space) onto
# the BANC EM voxel grid via:
#
#   Stage 0   LSM (Carl Zeiss multi-channel) -> NRRD (NC82 + GFP)
#             (inst/python/lsm_to_nrrd.py wrapping tifffile + SimpleITK)
#   Stage A'  native NC82 + GFP -> JRC2018U   (Elastix multi-resolution
#             affine + B-spline; param files in
#             inst/extdata/elastix_lm_to_jrc2018u/)
#   Stage B'' JRC2018U -> JRC2018F            (saalfeldlab
#             template-building RenderTransformed JAR + H5
#             nat.jrcbrains::JRC2018U_JRC2018F.h5; dfield is already
#             in the resampling direction we need, so no `-i`)
#   Stage C   JRC2018F -> BANC voxel space    (Elastix transformix +
#             bancr's bundled brain_240721 chain --- same as the
#             Kondo pipeline)
#   Stage D   BANC NRRD -> Neuroglancer precomputed
#
# Outputs a per-sample <gene>_<region>-<sex>.ng/ precomputed dir, a
# per-batch timings.csv + registry_entries.json, and emits the
# gsutil + jq commands needed to upload to the BANC LM bucket and
# merge into the master registry.
#
# Usage:
#   Rscript inst/scripts/deng_to_banc.R [test|brain|vnc|all]
#     test   --  4 CapaR samples (BRAIN-{F,M} + VNC-{F,M})
#     brain  -- every *BRAIN* file in the local mirror
#     vnc    -- every *VNC*   file in the local mirror (NB: VNC chain
#               needs jrc2018vnc -> BANC transformix params; not yet
#               wired up. Errors out for now.)
#     all    -- brain + vnc
#
# Resumable via existing-output detection (skips any
# <gene>_<sample>.ng/ already on disk).
#
# == Data source ==========================================================
#
# Deng, B., Li, Q., Liu, X., Cao, Y., Li, B., Qian, Y., Xu, R., Mao, R.,
# Zhou, E., Zhang, W., Huang, J., Rao, Y. (2019). Chemoconnectomics:
# Mapping Chemical Transmission in Drosophila. Neuron 101(5):876-893.e4.
# PMID 30709658. https://doi.org/10.1016/j.neuron.2019.01.045
#
# Raw imagery (T2A-Gal4 + UAS-FRT-Stop-FRT-mCD8::GFP knock-ins for
# receptors / neuropeptides / neurotransmitter-related GPCRs) was
# acquired in the Bowen / Wei-Chung Allen Lee lab as part of an
# in-prep BANC LM atlas; the local mirror is at
#   imaging-CCT-Bowen/<gene>-<...>-<BRAIN|VNC>-<sex>.lsm
#
# == Output =========================================================
#
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
#                                      deng_et_al_2019/
# Indexed by the master registry at
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
#                                      registry.json
# Browseable from R via `bancr::banc_lm_links` and
# `bancr:::banc_lm_volumes()`.

suppressMessages({
  library(neuronbridger)
  library(jsonlite)
})

# --- knobs ------------------------------------------------------------

LSM_DIR  <- "/Users/asbates/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/imaging-CCT-Bowen"
OUT_ROOT <- "~/deng_to_banc"
DATASET  <- "deng_et_al_2019"

mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "test"

# --- file-name parser -----------------------------------------------

# Canonical pattern: <gene>-<...details...>-(BRAIN|VNC)-<F|M>(-<trial>)?.lsm
# Some files lack the trailing -F/-M (older / unsexed); skip those for
# now. Some contain spaces / non-ASCII --- handled.
parse_lsm <- function(fname) {
  bn <- sub("\\.lsm$", "", basename(fname), ignore.case = TRUE)
  m <- regmatches(bn, regexec("^(.*?)[-_].*?(BRAIN|VNC)([-_][FM])?", bn,
                              ignore.case = TRUE))[[1]]
  if (!length(m) || length(m) < 3L) return(NULL)
  gene   <- m[[2]]
  region <- toupper(m[[3]])
  sex    <- if (length(m) >= 4L && nzchar(m[[4]])) toupper(sub("[-_]", "", m[[4]])) else "U"
  list(gene = gene, region = region, sex = sex,
       sample = sprintf("%s-%s", region, sex))
}

# --- target set per mode -------------------------------------------

GENE_SET_TEST <- c(
  "CapaR-L-T1-BRAIN-F", "CapaR-L-T1-BRAIN-M",
  "CapaR-L-T1-VNC-F",   "CapaR-L-T1-VNC-M"
)

list_lsms <- function(mode) {
  if (mode == "test") {
    sapply(GENE_SET_TEST, function(stem) file.path(LSM_DIR, paste0(stem, ".lsm")))
  } else {
    pat <- switch(mode,
                  brain = "BRAIN.*\\.lsm$",
                  vnc   = "VNC.*\\.lsm$",
                  all   = "(BRAIN|VNC).*\\.lsm$",
                  stop("mode must be test|brain|vnc|all; got: ", mode))
    list.files(LSM_DIR, pattern = pat, full.names = TRUE,
               recursive = FALSE, ignore.case = TRUE)
  }
}

if (!dir.exists(LSM_DIR))
  stop("Deng LSM dir not found: ", LSM_DIR)

out_root <- normalizePath(path.expand(OUT_ROOT), mustWork = FALSE)
out_dir  <- file.path(out_root, DATASET)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

targets <- list_lsms(mode[1])
cat(sprintf("== %d LSM target(s) for mode '%s' ==\n", length(targets), mode[1]))

results <- list()
timings <- list()
for (src in targets) {
  if (!file.exists(src)) {
    message(sprintf("SKIP  %s --- file not found", basename(src)))
    next
  }
  meta <- parse_lsm(src)
  if (is.null(meta)) {
    message(sprintf("SKIP  %s --- could not parse gene/region/sex", basename(src)))
    next
  }
  if (meta$region == "VNC") {
    message(sprintf("SKIP  %s --- VNC chain not yet wired (needs jrc2018vnc Stage C)",
                    basename(src)))
    next
  }
  gene <- meta$gene; sample <- meta$sample
  name <- paste(gene, sample, sep = "_")
  pc_dir <- file.path(out_dir, paste0(name, ".ng"))
  if (dir.exists(pc_dir) && length(list.files(pc_dir))) {
    message(sprintf("SKIP  %s --- precomputed dir already exists", name))
    next
  }
  message(sprintf("\n--- %s ---", name))
  message("  source: ", basename(src))

  res <- try(lsm_to_banc_layer(
    input        = src,
    gene         = gene,
    sample       = sample,
    channel      = "GFP",
    dataset      = DATASET,
    output_dir   = out_dir,
    source_path  = src,
    keep_intermediates = FALSE
  ), silent = TRUE)
  if (inherits(res, "try-error")) {
    message("  FAILED:\n  ", attr(res, "condition")$message)
    next
  }
  results[[name]] <- res
  timings[[length(timings) + 1L]] <- data.frame(
    name           = name,
    elastix_metric = res$elastix_metric,
    stage0         = res$timings[["stage0_lsm_extract"]],
    stageA         = res$timings[["stageA_native_to_jrcu"]],
    stageB         = res$timings[["stageB_jrcu_to_jrcf"]],
    stageC         = res$timings[["stageC_jrcf_to_banc"]],
    stageD         = res$timings[["stageD_precomputed"]],
    total          = res$timings[["total"]]
  )
}

if (!length(results)) {
  message("\nNo new volumes processed.")
  quit(save = "no", status = 0L)
}

times_df <- do.call(rbind, timings)
write.csv(times_df, file.path(out_dir, "timings.csv"), row.names = FALSE)
cat("\n=== timings (s) + Elastix final metric ===\n")
print(times_df, row.names = FALSE)

reg_entries <- do.call(rbind, lapply(results, `[[`, "registry_entry"))
reg_json    <- list(
  schema_version = 1L,
  volumes        = lapply(seq_len(nrow(reg_entries)), function(i) {
    row <- as.list(reg_entries[i, ])
    row$voxdims_nm <- as.integer(unlist(row$voxdims_nm))
    row
  })
)
writeLines(jsonlite::toJSON(reg_json, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "registry_entries.json"))

cat("\nNext steps (lee-lab maintainers only):\n")
cat(sprintf("  bash <upload helper> %s/*.ng %s\n", out_dir, DATASET))
cat("Then merge per-batch registry_entries.json into master registry:\n")
cat("  gsutil cat gs://.../light_level/registry.json | \\\n")
cat(sprintf("    jq '.volumes += $batch.volumes' --argfile batch %s/registry_entries.json | \\\n",
            out_dir))
cat("    gsutil cp - gs://.../light_level/registry.json\n")
