#!/usr/bin/env Rscript
#
# deng_brain_to_banc.R --- register every Deng et al. 2019 BRAIN LSM into
# BANC voxel space and upload the resulting Neuroglancer precomputed
# layers to GCS.
#
# Per-sample pipeline (orchestrated by neuronbridger::lsm_to_banc_layer):
#   Stage 0   LSM -> NC82 + GFP NRRDs                 (~5 s)
#   Stage A'  native -> JRC2018U via Elastix          (~90 s)
#   Stage B'' JRC2018U -> JRC2018F via H5             (~3 min)
#   Stage C   JRC2018F -> BANC voxel via transformix  (~3 min)
#   Stage D   BANC NRRD -> Neuroglancer precomputed   (~30 s)
#
# Outputs per sample <gene>_BRAIN-<sex>.ng/ at OUT_ROOT/<DATASET>/. After
# all samples finish, this script optionally uploads each .ng dir to GCS
# (gsutil cp -Z), merges the per-sample registry entries into the master
# registry.json, and re-mints bancr::banc_lm_links via the bancr
# data-raw/make_banc_lm_links.R helper.
#
# Usage:
#   Rscript inst/scripts/deng_brain_to_banc.R [test|all] [--upload] [--register]
#
#     test       --  4 CapaR test samples (BRAIN-{F,M} + VNC-{F,M};
#                    VNC samples are silently skipped here).
#     all        --  every *BRAIN*.lsm in the local mirror.
#     --upload   --  on completion, gsutil cp -Z the .ng dirs to
#                    gs://lee-lab.../light_level/<DATASET>/.
#     --register --  on completion, merge registry_entries.json into the
#                    master registry.json on GCS. Implies --upload.
#
# Re-runs are resumable: any sample whose .ng dir already exists locally
# is skipped.
#
# == Source =========================================================
# Deng, B. et al. (2019). Chemoconnectomics: Mapping Chemical Transmission
# in Drosophila. Neuron 101(5):876-893.e4. PMID 30709658.
# Raw imagery acquired in the Bowen / Lee lab; mirrored locally at
#   imaging-CCT-Bowen/<gene>-<...>-(BRAIN|VNC)-<sex>.lsm
#
# == Output =========================================================
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
#                                       deng_et_al_2019/<gene>_BRAIN-<sex>.ng/
# Indexed by:
#   gs://lee-lab.../light_level/registry.json
# Browseable via bancr::banc_lm_links / bancr:::banc_lm_volumes().

suppressMessages({
  library(neuronbridger)
  library(jsonlite)
})

# --- knobs ---------------------------------------------------------
LSM_DIR  <- "/Users/asbates/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/imaging-CCT-Bowen"
OUT_ROOT <- "~/deng_to_banc"
DATASET  <- "deng_et_al_2019"
GS_BASE  <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level"
BANCR_REPO <- "~/Projects/flyconnectome/bancr"

# --- args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("test", "all")) args[1] else "test"
upload   <- "--upload"   %in% args || "--register" %in% args
register <- "--register" %in% args

# --- file-name parser ----------------------------------------------
parse_lsm <- function(fname) {
  bn <- sub("\\.lsm$", "", basename(fname), ignore.case = TRUE)
  m <- regmatches(bn, regexec("^(.*?)[-_].*?(BRAIN|VNC)([-_][FM])?", bn,
                              ignore.case = TRUE))[[1]]
  if (!length(m) || length(m) < 3L) return(NULL)
  gene   <- m[[2]]
  region <- toupper(m[[3]])
  sex    <- if (length(m) >= 4L && nzchar(m[[4]]))
              toupper(sub("[-_]", "", m[[4]])) else "U"
  list(gene = gene, region = region, sex = sex,
       sample = sprintf("%s-%s", region, sex))
}

list_lsms <- function(mode) {
  if (mode == "test") {
    sapply(c("CapaR-L-T1-BRAIN-F", "CapaR-L-T1-BRAIN-M"),
           function(s) file.path(LSM_DIR, paste0(s, ".lsm")))
  } else {
    list.files(LSM_DIR, pattern = "BRAIN.*\\.lsm$", full.names = TRUE,
               recursive = FALSE, ignore.case = TRUE)
  }
}

# --- main loop -----------------------------------------------------
if (!dir.exists(LSM_DIR)) stop("LSM dir not found: ", LSM_DIR)
out_root <- normalizePath(path.expand(OUT_ROOT), mustWork = FALSE)
out_dir  <- file.path(out_root, DATASET)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

targets <- list_lsms(mode)
cat(sprintf("== %d BRAIN LSM target(s), mode '%s' ==\n", length(targets), mode))

results <- list(); timings <- list()
for (src in targets) {
  if (!file.exists(src)) {
    message("SKIP  ", basename(src), " --- file not found"); next
  }
  meta <- parse_lsm(src)
  if (is.null(meta) || meta$region != "BRAIN") next
  name   <- paste(meta$gene, meta$sample, sep = "_")
  pc_dir <- file.path(out_dir, paste0(name, ".ng"))
  if (dir.exists(pc_dir) && length(list.files(pc_dir))) {
    message("SKIP  ", name, " --- precomputed dir already exists"); next
  }
  message("\n--- ", name, " ---")
  res <- try(lsm_to_banc_layer(
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
write.csv(times_df, file.path(out_dir, "timings.csv"), row.names = FALSE)
cat("\n=== timings ===\n"); print(times_df, row.names = FALSE)

reg_entries <- do.call(rbind, lapply(results, `[[`, "registry_entry"))
reg_json <- list(
  schema_version = 1L,
  volumes = lapply(seq_len(nrow(reg_entries)), function(i) {
    r <- as.list(reg_entries[i, ])
    r$voxdims_nm <- as.integer(unlist(r$voxdims_nm)); r
  })
)
local_reg <- file.path(out_dir, "registry_entries.json")
writeLines(jsonlite::toJSON(reg_json, auto_unbox = TRUE, pretty = TRUE), local_reg)
cat("\nwrote", local_reg, "\n")

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
