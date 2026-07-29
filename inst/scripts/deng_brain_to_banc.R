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
#   Rscript inst/scripts/deng_brain_to_banc.R \
#     [test|all] [--upload] [--register] [--data-root <path>]
#
#     test            --  4 CapaR test samples (BRAIN-{F,M} + VNC-{F,M};
#                         VNC samples are silently skipped here).
#     all             --  every *BRAIN*.lsm in the local mirror.
#     --upload        --  on completion, gsutil cp -Z the .ng dirs to
#                         gs://lee-lab.../light_level/<DATASET>/.
#     --register      --  on completion, merge registry_entries.json into
#                         the master registry.json on GCS. Implies --upload.
#     --data-root P   --  override the local data root. Defaults to
#                         ~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat
#                         (or $NEURONBRIDGER_DATA_ROOT if set).
#     --one <name>    --  worker mode: process one sample (path or bare
#                         basename in LSM_DIR), write per-sample RDS
#                         stubs, exit. Used internally by outer mode so
#                         each sample runs in a fresh R heap (avoids
#                         mid-batch macOS memory-pressure slowdown).
#
# Local layout expected under <data-root>:
#   imaging-CCT-Bowen/   raw LSMs (sources)
#   templates/           JRC2018_UNISEX_20x_HR.nrrd, ...
#   deng_to_banc/        outputs (this script writes <DATASET>/ here)
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
# All local sources + outputs live under a single DATA_ROOT. The default
# is Alex's Dropbox layout; override via NEURONBRIDGER_DATA_ROOT or the
# --data-root CLI flag (handled below).
DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"
DATASET    <- "deng_et_al_2019"
GS_BASE    <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level"
BANCR_REPO <- "~/Projects/flyconnectome/bancr"

# --- args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
# pull --data-root <path> out of args if present
dr_idx <- which(args == "--data-root")
data_root <- if (length(dr_idx) && dr_idx < length(args)) args[dr_idx + 1L] else Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT)
DATA_ROOT <- path.expand(data_root)
if (length(dr_idx)) args <- args[-c(dr_idx, dr_idx + 1L)]

force    <- "--force" %in% args
upload   <- "--upload"   %in% args || "--register" %in% args
register <- "--register" %in% args
one_idx  <- which(args == "--one")
one      <- if (length(one_idx) && one_idx < length(args)) args[one_idx + 1L] else NA_character_
if (length(one_idx)) args <- args[-c(one_idx, one_idx + 1L)]
args     <- args[!args %in% c("--force", "--upload", "--register")]
mode     <- if (length(args) && args[1] %in% c("test", "all")) args[1] else "test"

LSM_DIR <- file.path(DATA_ROOT, "imaging-CCT-Bowen")
OUT_DIR <- file.path(DATA_ROOT, "deng_to_banc", DATASET)
TPL_DIR <- file.path(DATA_ROOT, "templates")
options(neuronbridger.jrc2018u_template =
          file.path(TPL_DIR, "JRC2018_UNISEX_20x_HR.nrrd"))
cat("DATA_ROOT: ", DATA_ROOT, "\n", sep = "")

# --- file-name parser ----------------------------------------------
parse_lsm <- function(fname) {
  bn <- sub("\\.lsm$", "", basename(fname), ignore.case = TRUE)
  # If the filename has multiple delimited BRAIN|VNC tokens (e.g.
  # "brain+vnc-vnc.lsm"), splice out the middle so the region tag we pick
  # up below is the trailing one. Word-bounded match so substrings inside
  # tokens like "VNCGFP" don't count.
  tag_pat <- "(^|[^A-Za-z])([Bb][Rr][Aa][Ii][Nn]|[Vv][Nn][Cc])(?![A-Za-z])"
  tag_matches <- gregexpr(tag_pat, bn, perl = TRUE)[[1]]
  if (length(tag_matches) > 1L && !any(tag_matches == -1L)) {
    first <- tag_matches[1L]; last <- tag_matches[length(tag_matches)]
    first_actual <- if (first == 1L) 1L else first + 1L
    last_actual  <- if (last  == 1L) 1L else last  + 1L
    bn <- paste0(substr(bn, 1L, first_actual - 2L), "-",
                 substr(bn, last_actual, nchar(bn)))
  }
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
out_dir <- normalizePath(OUT_DIR, mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stubs_dir <- file.path(out_dir, "_per_sample")
dir.create(stubs_dir, showWarnings = FALSE, recursive = TRUE)

# --- worker mode: one sample, fresh R heap, drop RDS stubs, exit ---
if (!is.na(one)) {
  src <- if (file.exists(one)) one else file.path(LSM_DIR, sub("\\.lsm$", "", basename(one), ignore.case = TRUE))
  if (!file.exists(src) && !grepl("\\.lsm$", src, ignore.case = TRUE))
    src <- paste0(src, ".lsm")
  if (!file.exists(src)) stop("--one target not found: ", one)
  meta <- parse_lsm(src)
  if (is.null(meta) || meta$region != "BRAIN") stop("--one target didn't parse as BRAIN: ", src)
  name <- paste(meta$gene, meta$sample, sep = "_")
  message("--- ", name, " (worker, PID=", Sys.getpid(), ") ---")
  res <- lsm_to_banc_layer(
    input              = src,
    gene               = meta$gene,
    sample             = meta$sample,
    channel            = "GFP",
    dataset            = DATASET,
    output_dir         = out_dir,
    source_path        = src,
    keep_intermediates = FALSE
  )
  saveRDS(data.frame(name = name,
                     elastix_metric = res$elastix_metric,
                     total          = res$timings[["total"]]),
          file.path(stubs_dir, paste0(name, ".timings.rds")))
  saveRDS(res$registry_entry,
          file.path(stubs_dir, paste0(name, ".registry.rds")))
  quit(save = "no", status = 0L)
}

targets <- list_lsms(mode)
cat(sprintf("== %d BRAIN LSM target(s), mode '%s' ==\n", length(targets), mode))

processed_names <- character(0); failed_names <- character(0)
self <- sub("^--file=", "",
            commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1L])
self <- normalizePath(self, mustWork = TRUE)

for (src in targets) {
  if (!file.exists(src)) {
    message("SKIP  ", basename(src), " --- file not found"); next
  }
  meta <- parse_lsm(src)
  if (is.null(meta) || meta$region != "BRAIN") next
  name   <- paste(meta$gene, meta$sample, sep = "_")
  pc_dir <- file.path(out_dir, paste0(name, ".ng"))
  if (dir.exists(pc_dir) && length(list.files(pc_dir)) && !force) {
    message("SKIP  ", name, " (precomputed dir exists; pass --force to recompute)")
    if (file.exists(file.path(stubs_dir, paste0(name, ".registry.rds"))))
      processed_names <- c(processed_names, name)
    next
  }
  worker_args <- c(shQuote(self), "--one", shQuote(src),
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

if (!length(results_reg)) { cat("\nNo new volumes processed.\n"); quit(save = "no") }
if (!length(timings)) { cat("\nNo timings — RDS stubs missing?\n"); quit(save = "no") }

times_df <- do.call(rbind, timings)
write.csv(times_df, file.path(out_dir, "timings.csv"), row.names = FALSE)
cat("\n=== timings ===\n"); print(times_df, row.names = FALSE)

reg_entries <- do.call(rbind, results_reg)
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
