#!/usr/bin/env Rscript
#
# deng_vnc_to_banc.R --- register every Deng et al. 2019 VNC LSM into BANC
# voxel space and upload the resulting Neuroglancer precomputed layers
# to GCS.
#
# Per-sample pipeline (orchestrated by neuronbridger::lsm_vnc_to_banc_layer):
#   Stage 0   LSM -> NC82 + GFP NRRDs                       (~5 s)
#   Stage A'  native -> JRC2018VNCF via Elastix             (~90 s)
#   Stage B   JRC2018VNCF -> BANC voxel via tpsreg          (~30 min,
#             inverse-warp + max-pool resample; uses
#             bancr::banc_to_jrcvnc2018f_tpsreg via the Python helper
#             at inst/python/lm_vnc_inverse_warp.py).
#   Stage D   BANC NRRD -> Neuroglancer precomputed         (~30 s)
#
# The bancr `vnc_240721/BANC_to_template.txt` Elastix bridge is
# incomplete (only inverts the fine bspline; the manual + elastix
# affines are missing) so we go via the bundled tpsreg instead. The
# tpsreg is mathematically equivalent to bancr::banc_to_JRC2018F(...,
# region = "vnc", inverse = TRUE) and lands signal at the BANC VNC
# neuropil mesh's coordinate frame.
#
# Outputs per sample <gene>_VNC-<sex>.ng/ at OUT_ROOT/<DATASET>/.
# After all samples finish, this script optionally uploads each .ng
# dir to GCS, merges per-sample registry entries into the master
# registry.json, and re-mints bancr::banc_lm_links.
#
# Usage:
#   Rscript inst/scripts/deng_vnc_to_banc.R \
#     [test|all] [--upload] [--register] [--data-root <path>]
#
#     test            --  CapaR-L-T1-VNC-{F,M}
#     all             --  every *VNC*.lsm in the local mirror
#     --upload        --  on completion, gsutil cp -Z .ng dirs to GCS
#     --register      --  on completion, merge registry + remint banc_lm_links.
#                         Implies --upload.
#     --data-root P   --  override the local data root. Defaults to
#                         ~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat
#                         (or $NEURONBRIDGER_DATA_ROOT if set).
#
# Local layout expected under <data-root>:
#   imaging-CCT-Bowen/   raw LSMs (sources)
#   templates/           JRC2018_VNC_FEMALE_461.nrrd, ...
#   deng_to_banc/        outputs (this script writes <DATASET>/ here)
#
# Re-runs are resumable: any sample whose .ng dir already exists is skipped.

suppressMessages({
  library(neuronbridger)
  library(jsonlite)
})

# --- knobs ---------------------------------------------------------
DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"
DATASET    <- "deng_et_al_2019"
GS_BASE    <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level"
BANCR_REPO <- "~/Projects/flyconnectome/bancr"

# --- args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
dr_idx <- which(args == "--data-root")
data_root <- if (length(dr_idx) && dr_idx < length(args)) args[dr_idx + 1L] else Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT)
DATA_ROOT <- path.expand(data_root)
if (length(dr_idx)) args <- args[-c(dr_idx, dr_idx + 1L)]

force    <- "--force" %in% args
upload   <- "--upload"   %in% args || "--register" %in% args
register <- "--register" %in% args
args     <- setdiff(args, c("--force", "--upload", "--register"))
mode     <- if (length(args) && args[1] %in% c("test", "all")) args[1] else "test"

LSM_DIR <- file.path(DATA_ROOT, "imaging-CCT-Bowen")
OUT_DIR <- file.path(DATA_ROOT, "deng_to_banc", DATASET)
TPL_DIR <- file.path(DATA_ROOT, "templates")
options(neuronbridger.jrcvnc2018f_template =
          file.path(TPL_DIR, "JRC2018_VNC_FEMALE_461.nrrd"))
cat("DATA_ROOT: ", DATA_ROOT, "\n", sep = "")

parse_lsm <- function(fname) {
  bn <- sub("\\.lsm$", "", basename(fname), ignore.case = TRUE)
  m <- regmatches(bn, regexec("^(.*?)[-_].*?(BRAIN|VNC)([-_][FM])?", bn,
                              ignore.case = TRUE))[[1]]
  if (!length(m) || length(m) < 3L) return(NULL)
  gene <- m[[2]]; region <- toupper(m[[3]])
  sex  <- if (length(m) >= 4L && nzchar(m[[4]]))
            toupper(sub("[-_]", "", m[[4]])) else "U"
  list(gene = gene, region = region, sex = sex,
       sample = sprintf("%s-%s", region, sex))
}

list_lsms <- function(mode) {
  if (mode == "test") {
    sapply(c("CapaR-L-T1-VNC-F", "CapaR-L-T1-VNC-M"),
           function(s) file.path(LSM_DIR, paste0(s, ".lsm")))
  } else {
    list.files(LSM_DIR, pattern = "VNC.*\\.lsm$", full.names = TRUE,
               recursive = FALSE, ignore.case = TRUE)
  }
}

# --- main ----------------------------------------------------------
if (!dir.exists(LSM_DIR)) stop("LSM dir not found: ", LSM_DIR)
out_dir <- normalizePath(OUT_DIR, mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

targets <- list_lsms(mode)
cat(sprintf("== %d VNC LSM target(s), mode '%s' ==\n", length(targets), mode))

results <- list(); timings <- list()
for (src in targets) {
  if (!file.exists(src)) { message("SKIP  ", basename(src), " (not found)"); next }
  meta <- parse_lsm(src)
  if (is.null(meta) || meta$region != "VNC") next
  name   <- paste(meta$gene, meta$sample, sep = "_")
  pc_dir <- file.path(out_dir, paste0(name, ".ng"))
  if (dir.exists(pc_dir) && length(list.files(pc_dir)) && !force) {
    message("SKIP  ", name, " (precomputed dir exists; pass --force to recompute)"); next
  }
  message("\n--- ", name, " ---")
  res <- try(lsm_vnc_to_banc_layer(
    input              = src,
    gene               = meta$gene,
    sample             = meta$sample,
    channel            = "GFP",
    dataset            = DATASET,
    output_dir         = out_dir,
    source_path        = src,
    keep_intermediates = FALSE
  ), silent = TRUE)
  if (inherits(res, "try-error")) {
    message("  FAILED: ", attr(res, "condition")$message); next
  }
  results[[name]] <- res
  timings[[length(timings) + 1L]] <- data.frame(
    name = name, elastix_metric = res$elastix_metric, total = res$timings[["total"]]
  )
}

if (!length(results)) { cat("\nNo new volumes processed.\n"); quit(save = "no") }

times_df <- do.call(rbind, timings)
write.csv(times_df, file.path(out_dir, "timings_vnc.csv"), row.names = FALSE)
cat("\n=== timings ===\n"); print(times_df, row.names = FALSE)

reg_entries <- do.call(rbind, lapply(results, `[[`, "registry_entry"))
reg_json <- list(
  schema_version = 1L,
  volumes = lapply(seq_len(nrow(reg_entries)), function(i) {
    r <- as.list(reg_entries[i, ])
    r$voxdims_nm <- as.integer(unlist(r$voxdims_nm)); r
  })
)
local_reg <- file.path(out_dir, "registry_entries_vnc.json")
writeLines(jsonlite::toJSON(reg_json, auto_unbox = TRUE, pretty = TRUE), local_reg)
cat("\nwrote", local_reg, "\n")

if (upload) {
  cat("\n=== uploading .ng dirs to GCS ===\n")
  ng_paths <- file.path(out_dir, paste0(names(results), ".ng"))
  ng_paths <- ng_paths[dir.exists(ng_paths)]
  cmd <- sprintf("gsutil -m cp -Z -r %s %s/%s/",
                 paste(shQuote(ng_paths), collapse = " "),
                 GS_BASE, DATASET)
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
