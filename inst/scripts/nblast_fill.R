#!/usr/bin/env Rscript
#
# nblast_fill.R --- fill in NBLAST + voxel_attr columns of existing
# `_hits.csv` files produced by colormip_search_against_banc.R. Skips
# rows that already have both scores populated (safe to re-run).
#
# The colormip_search_against_banc.R pass sets NBLAST and voxel_attr
# to NA whenever the top-K hit's L2 SWC is not in the local mirror
# under <data-root>/banc_meta/banc_banc_space_swc/. This helper:
#   1. reads every _hits.csv under HITS_DIR
#   2. gathers the union of preferred_ids
#   3. bulk-`gsutil cp` any missing <id>_l2.swc from GCS
#   4. re-computes NBLAST + voxel_attr per (query, hit) using the
#      SAME code path colormip_search_against_banc.R uses
#   5. rewrites each _hits.csv in place and refreshes lm_to_banc_colormip_hits.csv
#
# Usage:
#   Rscript inst/scripts/nblast_fill.R REGION HITS_DIR WARPS_ROOT [--data-root P]
#
#     REGION       "brain" | "VNC"
#     HITS_DIR     dir of <sample>_hits.csv from an earlier search run
#     WARPS_ROOT   dir containing <sample>.nrrd (same as make_colormip_library
#                  output_dir for the region)
#     --data-root  local data root (default: NEURONBRIDGER_DATA_ROOT or Dropbox)

suppressMessages({
  pkg_dir <- normalizePath(getwd(), mustWork = FALSE)
  if (file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      identical(read.dcf(file.path(pkg_dir, "DESCRIPTION"), "Package")[[1]],
                "neuronbridger")) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else library(neuronbridger)
  library(bancr); library(nat); library(nat.flybrains)
  library(nat.templatebrains); library(nat.jrcbrains); library(nat.nblast)
  nat.jrcbrains::register_saalfeldlab_registrations()
  library(readr); library(dplyr); library(reticulate); library(FNN)
})

# --- args ---------------------------------------------------------
DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"
args <- commandArgs(trailingOnly = TRUE)
dr_i <- which(args == "--data-root")
data_root <- if (length(dr_i) && dr_i < length(args)) {
  args[dr_i + 1L]
} else {
  Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT)
}
if (length(dr_i)) args <- args[-c(dr_i, dr_i + 1L)]

if (length(args) < 3) stop("usage: nblast_fill.R REGION HITS_DIR WARPS_ROOT [--data-root P]")
REGION     <- match.arg(args[1], c("brain", "VNC"))
HITS_DIR   <- path.expand(args[2])
WARPS_ROOT <- path.expand(args[3])
DATA_ROOT  <- path.expand(data_root)
SWC_DIR    <- file.path(DATA_ROOT, "banc_meta", "banc_banc_space_swc")
GS_SWC     <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_banc_space_swc"
dir.create(SWC_DIR, showWarnings = FALSE, recursive = TRUE)

# --- phase 1: gather all preferred_ids across hit CSVs -----------
csvs <- list.files(HITS_DIR, pattern = "_hits\\.csv$", full.names = TRUE)
cat(sprintf("== reading %d hit CSVs ==\n", length(csvs)))
all_ids <- unique(unlist(lapply(csvs, function(f) {
  d <- readr::read_csv(f, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  d$preferred_id
})))
cat(sprintf("  union preferred_ids: %d unique\n", length(all_ids)))

# --- phase 2: bulk-sync missing SWCs -----------------------------
have_local <- basename(list.files(SWC_DIR, "_l2\\.swc$"))
have_local <- sub("_l2\\.swc$", "", have_local)
missing <- setdiff(all_ids, have_local)
cat(sprintf("  local SWCs: %d ; missing: %d\n", length(have_local), length(missing)))
if (length(missing)) {
  tmp_list <- tempfile(fileext = ".txt")
  writeLines(sprintf("%s/%s_l2.swc", GS_SWC, missing), tmp_list)
  cmd <- sprintf("cat %s | gsutil -m cp -I %s/", shQuote(tmp_list), shQuote(SWC_DIR))
  cat("  running: ", cmd, "\n"); t0 <- Sys.time()
  rc <- system(cmd)
  cat(sprintf("  bulk-cp exit=%d, %.1f min\n", rc,
              as.numeric(difftime(Sys.time(), t0, units="mins"))))
  have_local <- basename(list.files(SWC_DIR, "_l2\\.swc$"))
  have_local <- sub("_l2\\.swc$", "", have_local)
  cat(sprintf("  local SWCs after sync: %d\n", length(have_local)))
}

# --- python helper (trilinear sample of warped LM) ---------------
PYBIN <- Sys.getenv("RETICULATE_PYTHON",
                    unset = "/Users/asbates/Library/r-miniconda-arm64/envs/r-reticulate/bin/python3")
reticulate::use_python(PYBIN, required = TRUE)
sitk <- reticulate::import("SimpleITK", convert = FALSE)
np   <- reticulate::import("numpy",     convert = FALSE)

# --- phase 3: per-CSV fill in -----------------------------------
process_one <- function(csv_path) {
  name <- sub("_hits\\.csv$", "", basename(csv_path))
  d <- readr::read_csv(csv_path, show_col_types = FALSE,
                       col_types = readr::cols(banc_id = "c", root_888 = "c",
                                                preferred_id = "c"))
  # short-circuit if all NBLAST + voxel_attr already populated
  need_nb <- is.na(d$nblast)
  need_va <- is.na(d$voxel_attr)
  if (!any(need_nb | need_va)) {
    cat(sprintf("SKIP %s (all NBLAST + voxel_attr already populated)\n", name))
    return(invisible(NULL))
  }
  warp_nrrd <- file.path(WARPS_ROOT, paste0(name, ".nrrd"))
  if (!file.exists(warp_nrrd)) {
    cat(sprintf("WARN %s: no warped LM at %s -- skipping\n", name, warp_nrrd))
    return(invisible(NULL))
  }
  cat(sprintf("\n=== %s (%d rows, %d need NBLAST, %d need voxel_attr) ===\n",
              name, nrow(d), sum(need_nb), sum(need_va)))

  img <- sitk$ReadImage(warp_nrrd)
  arr <- np$clip(sitk$GetArrayFromImage(img), 0L, 255L)$astype(np$uint8)
  sp  <- unlist(reticulate::py_to_r(img$GetSpacing()))
  zyx <- np$nonzero(np$greater(arr, np$uint8(20)))
  zs  <- reticulate::py_to_r(np$asarray(zyx[[0]]))
  ys  <- reticulate::py_to_r(np$asarray(zyx[[1]]))
  xs  <- reticulate::py_to_r(np$asarray(zyx[[2]]))
  n   <- length(xs)
  if (!n) { cat("  no LM foreground; leaving as NA\n"); return(invisible(NULL)) }
  set.seed(1); idx <- if (n > 50000L) sample.int(n, 50000L) else seq_len(n)
  lm_pts_um <- cbind(xs[idx]*sp[1], ys[idx]*sp[2], zs[idx]*sp[3])
  lm_dp <- nat::dotprops(lm_pts_um, k = 5L)

  target_tb  <- if (REGION == "brain") "JRC2018U"  else "JRCVNC2018U"
  source_tb  <- if (REGION == "brain") "JRC2018F"  else "JRCVNC2018F"
  region_tps <- if (REGION == "brain") "brain"     else "vnc"

  # NBLAST per hit
  d$nblast <- vapply(seq_len(nrow(d)), function(i) {
    swc <- file.path(SWC_DIR, paste0(d$preferred_id[i], "_l2.swc"))
    if (!file.exists(swc)) return(NA_real_)
    tryCatch({
      sk <- nat::read.neuron(swc)
      pts <- nat::xyzmatrix(sk); if (nrow(pts) < 6L) return(NA_real_)
      pjf <- bancr::banc_to_JRC2018F(pts, method = "tpsreg",
                                     banc.units = "nm", region = region_tps)
      ptg <- nat.templatebrains::xform_brain(pjf, sample = source_tb,
                                             reference = target_tb)
      if (nrow(ptg) < 6L) return(NA_real_)
      dp <- nat::dotprops(ptg, k = 5L)
      as.numeric(nat.nblast::nblast(dp, target = lm_dp, normalised = TRUE))
    }, error = function(e) NA_real_)
  }, numeric(1))

  # voxel_attr: pool SWCs → NN
  tps_fwd <- if (REGION == "VNC") bancr::jrcvnc2018f_to_banc_tpsreg
             else                  bancr::jrc2018f_to_banc_tpsreg
  lm_nm <- Morpho::tps3d(lm_pts_um, refmat = tps_fwd$refmat,
                         tarmat = tps_fwd$tarmat, lambda = tps_fwd$lambda,
                         threads = 4)
  ptlist <- list(); lab <- list()
  for (i in seq_len(nrow(d))) {
    swc <- file.path(SWC_DIR, paste0(d$preferred_id[i], "_l2.swc"))
    if (!file.exists(swc)) next
    sk <- tryCatch(nat::read.neuron(swc), error = function(e) NULL)
    if (is.null(sk)) next
    p <- nat::xyzmatrix(sk); if (!nrow(p)) next
    ptlist[[length(ptlist)+1]] <- p; lab[[length(lab)+1]] <- rep(i, nrow(p))
  }
  if (length(ptlist)) {
    P <- do.call(rbind, ptlist); L <- unlist(lab)
    nn <- FNN::get.knnx(P, lm_nm, k = 1)
    votes <- L[nn$nn.index[,1]]
    d$voxel_attr <- tabulate(votes, nbins = nrow(d)) / length(votes)
  }

  # rerank
  d$cm_rank <- dplyr::min_rank(dplyr::desc(d$cm_score))
  d$nb_rank <- dplyr::min_rank(dplyr::desc(replace(d$nblast, is.na(d$nblast), -Inf)))
  d$va_rank <- dplyr::min_rank(dplyr::desc(replace(d$voxel_attr, is.na(d$voxel_attr), -Inf)))
  d$mean_rank <- (d$cm_rank + d$nb_rank + d$va_rank) / 3

  d <- dplyr::arrange(d, mean_rank)
  readr::write_csv(d, csv_path)
  cat(sprintf("  wrote %s (%d NBLAST filled, %d voxel_attr>0)\n",
              csv_path, sum(!is.na(d$nblast)),
              sum(!is.na(d$voxel_attr) & d$voxel_attr > 0)))
}

for (f in csvs) process_one(f)

# --- refresh lm_to_banc_colormip_hits.csv ------------------------
csvs <- list.files(HITS_DIR, pattern = "_hits\\.csv$", full.names = TRUE)
master <- dplyr::bind_rows(lapply(csvs, function(f)
  readr::read_csv(f, show_col_types = FALSE,
                  col_types = readr::cols(banc_id = "c", root_888 = "c",
                                          preferred_id = "c"))))
master_path <- file.path(HITS_DIR, "lm_to_banc_colormip_hits.csv")
readr::write_csv(master, master_path)
cat(sprintf("\nwrote %d rows -> %s\n", nrow(master), master_path))
