#!/usr/bin/env Rscript
#
# colormip_search_against_banc.R --- search a directory of query
# colour-MIPs against the BANC neuron MIP library and rank hits by
# colour-MIP overlap, NBLAST morphological similarity, and per-voxel
# attribution to the closest BANC L2 skeleton.
#
# Pipeline per query:
#   1. colormip_search() vs the BANC library at
#      gs://lee-lab.../neuron_colormips/template_alignment_240721/
#      (locally synced to ~/banc_colormips/<template>/) and keep top-K.
#   2. Look up root_888 + cell metadata from
#      ~/banc_meta/banc_888_meta.feather.
#   3. NBLAST top-K BANC neurons against the LM-derived dotprops
#      using locally cached L2 SWCs at
#      ~/banc_meta/banc_banc_space_swc/<root_888>_l2.swc (no CAVE
#      backend hits; pre-fetch SWCs once for unique IDs across all
#      samples — see prepare_skeletons() below).
#   4. Voxel-attribution: forward-warp each LM voxel via tpsreg into
#      BANC nm and assign it to its nearest L2 skeleton; per-hit
#      vote-fraction is the third score column.
#   5. Combine cm_score, nblast, voxel_attr into mean_rank.
#
# Output: a single master CSV at <OUT>/master_augmented.csv with one
# row per (query, BANC hit) carrying all three scores plus
# super_cluster / cell_function / neurotransmitter_predicted.
#
# Usage:
#   Rscript inst/scripts/colormip_search_against_banc.R \
#     <query_mip_dir> <region (brain|VNC)> <warps_root> <out_dir> \
#     [TOP_K] [--force] [--data-root <path>]
#
#     query_mip_dir   directory of query colour-MIP PNGs (one per
#                     sample); typically the output of
#                     make_colormip_library.R.
#     region          "brain" -> JRC2018_UNISEX_20x_HR colormip lib;
#                     "VNC"   -> JRC2018_VNC_UNISEX_461 colormip lib.
#     warps_root      directory containing the JRC2018U_HR /
#                     JRCVNC2018U_HR-aligned LM volumes
#                     (<sample>_in_<template>.nrrd) — needed for the
#                     LM dotprops + voxel-attribution steps. Pass the
#                     "elastix" subdir or wherever your warped GFP
#                     NRRDs live; the script looks for
#                     <warps_root>/<sample>/elastix/gfp_xform/result.nrrd
#                     and <warps_root>/<sample>.nrrd.
#     out_dir         per-sample <name>_hits.csv files + master CSV
#                     are written here.
#     TOP_K           per-sample top hits to keep (default 25).
#     --force         recompute even if <name>_hits.csv already exists.
#     --data-root P   override local data root (used to find the BANC
#                     colorMIP library + v888 metadata + L2 SWCs).
#                     Default: ~/Library/CloudStorage/Dropbox-HMS/...
#
# Local layout expected under <data-root>:
#   banc_colormips/{JRC2018_UNISEX_20x_HR,JRC2018_VNC_UNISEX_461}/
#                            BANC neuron colorMIP libraries
#   banc_meta/banc_888_meta.feather             v888 metadata
#   banc_meta/banc_banc_space_swc/<id>_l2.swc   v888 L2 skeletons

suppressMessages({
  pkg_dir <- normalizePath(getwd(), mustWork = FALSE)
  if (file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      identical(read.dcf(file.path(pkg_dir, "DESCRIPTION"), "Package")[[1]],
                "neuronbridger")) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    library(neuronbridger)
  }
  library(bancr); library(arrow); library(nat); library(nat.flybrains)
  library(nat.templatebrains); library(nat.jrcbrains); library(nat.nblast)
  nat.jrcbrains::register_saalfeldlab_registrations()
  library(readr); library(dplyr); library(tibble); library(reticulate)
  library(Morpho); library(FNN)
})

# --- args ----------------------------------------------------------
DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"
args <- commandArgs(trailingOnly = TRUE)
dr_idx <- which(args == "--data-root")
data_root <- if (length(dr_idx) && dr_idx < length(args)) args[dr_idx + 1L] else Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT)
DATA_ROOT <- path.expand(data_root)
if (length(dr_idx)) args <- args[-c(dr_idx, dr_idx + 1L)]
force <- "--force" %in% args
args  <- args[args != "--force"]

if (length(args) < 4) stop("usage: colormip_search_against_banc.R QUERY_DIR REGION WARPS_ROOT OUT_DIR [TOP_K] [--force] [--data-root <path>]")
QUERY_DIR  <- path.expand(args[1])
REGION     <- match.arg(args[2], c("brain", "VNC"))
WARPS_ROOT <- path.expand(args[3])
OUT_DIR    <- path.expand(args[4])
TOP_K      <- if (length(args) >= 5) as.integer(args[5]) else 25L
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("DATA_ROOT: ", DATA_ROOT, "\n", sep = "")

LIB_DIR <- file.path(DATA_ROOT, "banc_colormips",
                     if (REGION == "brain") "JRC2018_UNISEX_20x_HR"
                     else "JRC2018_VNC_UNISEX_461")
META    <- file.path(DATA_ROOT, "banc_meta", "banc_888_meta.feather")
SWC_DIR <- file.path(DATA_ROOT, "banc_meta", "banc_banc_space_swc")
GS_LIB  <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_colormips/template_alignment_240721"
GS_SWC  <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_banc_space_swc"

if (!dir.exists(LIB_DIR) || !length(list.files(LIB_DIR, "\\.png$")))
  stop("BANC ", REGION, " MIP library missing at ", LIB_DIR,
       ". Run: gsutil -m rsync -r ", GS_LIB, "/",
       if (REGION == "brain") "JRC2018_UNISEX_20x_HR/" else "JRC2018_VNC_UNISEX_461/",
       " ", LIB_DIR, "/")

# --- python (for trilinear sample of the warped LM) ----------------
PYBIN <- Sys.getenv("RETICULATE_PYTHON",
                    unset = "/Users/asbates/Library/r-miniconda-arm64/envs/r-reticulate/bin/python3")
reticulate::use_python(PYBIN, required = TRUE)
sitk <- reticulate::import("SimpleITK", convert = FALSE)
np   <- reticulate::import("numpy",     convert = FALSE)

# --- v888 metadata -------------------------------------------------
m888 <- arrow::read_feather(path.expand(META))
m888$root_626 <- as.character(m888$root_626)
m888$root_888 <- as.character(m888$root_888)
META_LK <- m888 |> dplyr::select(root_626, root_888, super_cluster,
                                 cell_function, neurotransmitter_predicted) |>
  dplyr::distinct(root_626, .keep_all = TRUE)

# --- per-query worker ----------------------------------------------
process_query <- function(mip_path) {
  name <- sub("\\.png$", "", basename(mip_path))
  out_csv <- file.path(OUT_DIR, paste0(name, "_hits.csv"))
  if (file.exists(out_csv) && !force) {
    message("SKIP ", name, " (CSV exists; pass --force to recompute)")
    return(out_csv)
  }
  message("\n=== ", name, " ===")

  # 1. colorMIP search
  hits <- colormip_search(query = mip_path, library = LIB_DIR,
                          threshold = 100L, z_tolerance = 2L, xy_shift = 2L,
                          mirror = TRUE, top_n = TOP_K, mc.cores = 8L,
                          verbose = FALSE)
  if (!nrow(hits)) { message("no hits"); return(NULL) }
  hit_ids <- as.character(sub(".*/(\\d+)_in_.*", "\\1", hits$path))

  # 2. v888 lookup
  hits$banc_id <- hit_ids
  aug <- dplyr::left_join(hits, META_LK, by = c("banc_id" = "root_626"))
  aug$preferred_id <- ifelse(is.na(aug$root_888), aug$banc_id, aug$root_888)

  # 3. + 4. need the warped LM volume.
  warp_nrrd <- file.path(WARPS_ROOT, paste0(name, ".nrrd"))
  if (!file.exists(warp_nrrd))
    warp_nrrd <- file.path(WARPS_ROOT, paste0(name, "_in_",
                                              if (REGION == "brain") "JRC2018U_HR"
                                              else "JRCVNC2018U_HR", ".nrrd"))
  if (!file.exists(warp_nrrd)) {
    message("WARN: no warped LM at ", warp_nrrd, " — NBLAST + voxel_attr stay NA")
    aug$nblast <- NA_real_; aug$voxel_attr <- NA_real_
  } else {
    img <- sitk$ReadImage(warp_nrrd)
    arr <- np$clip(sitk$GetArrayFromImage(img), 0L, 255L)$astype(np$uint8)
    sp  <- unlist(reticulate::py_to_r(img$GetSpacing()))
    zyx <- np$nonzero(np$greater(arr, np$uint8(20)))
    zs  <- reticulate::py_to_r(np$asarray(zyx[[0]]))
    ys  <- reticulate::py_to_r(np$asarray(zyx[[1]]))
    xs  <- reticulate::py_to_r(np$asarray(zyx[[2]]))
    n <- length(xs)
    if (n) {
      set.seed(1); idx <- if (n > 50000L) sample.int(n, 50000L) else seq_len(n)
      lm_pts_um <- cbind(xs[idx]*sp[1], ys[idx]*sp[2], zs[idx]*sp[3])
      lm_dp <- nat::dotprops(lm_pts_um, k = 5L)

      # NBLAST per-hit using local SWC (no CAVE backend)
      target_tb <- if (REGION == "brain") "JRC2018U"  else "JRCVNC2018U"
      source_tb <- if (REGION == "brain") "JRC2018F"  else "JRCVNC2018F"
      region_tps <- if (REGION == "brain") "brain" else "vnc"

      aug$nblast <- vapply(seq_len(nrow(aug)), function(i) {
        swc <- file.path(SWC_DIR, paste0(aug$preferred_id[i], "_l2.swc"))
        if (!file.exists(swc)) return(NA_real_)
        tryCatch({
          sk <- nat::read.neuron(swc)
          pts <- nat::xyzmatrix(sk)
          if (nrow(pts) < 6L) return(NA_real_)
          pjf <- bancr::banc_to_JRC2018F(pts, method = "tpsreg",
                                         banc.units = "nm", region = region_tps)
          ptg <- nat.templatebrains::xform_brain(pjf, sample = source_tb,
                                                 reference = target_tb)
          if (nrow(ptg) < 6L) return(NA_real_)
          dp <- nat::dotprops(ptg, k = 5L)
          as.numeric(nat.nblast::nblast(dp, target = lm_dp, normalised = TRUE))
        }, error = function(e) NA_real_)
      }, numeric(1))

      # voxel-attribution: forward-warp LM -> BANC nm, NN to top-K SWCs
      tps_fwd <- if (REGION == "VNC") bancr::jrcvnc2018f_to_banc_tpsreg
                 else                  bancr::jrc2018f_to_banc_tpsreg
      lm_nm <- Morpho::tps3d(lm_pts_um, refmat = tps_fwd$refmat,
                             tarmat = tps_fwd$tarmat, lambda = tps_fwd$lambda,
                             threads = 4)
      ptlist <- list(); lab <- list()
      for (i in seq_len(nrow(aug))) {
        swc <- file.path(SWC_DIR, paste0(aug$preferred_id[i], "_l2.swc"))
        if (!file.exists(swc)) next
        sk <- tryCatch(nat::read.neuron(swc), error = function(e) NULL)
        if (is.null(sk)) next
        p <- nat::xyzmatrix(sk); if (!nrow(p)) next
        ptlist[[length(ptlist)+1]] <- p
        lab[[length(lab)+1]]       <- rep(i, nrow(p))
      }
      if (length(ptlist)) {
        P <- do.call(rbind, ptlist); L <- unlist(lab)
        nn <- FNN::get.knnx(P, lm_nm, k = 1)
        votes <- L[nn$nn.index[,1]]
        aug$voxel_attr <- tabulate(votes, nbins = nrow(aug)) / length(votes)
      } else aug$voxel_attr <- NA_real_
    } else { aug$nblast <- NA_real_; aug$voxel_attr <- NA_real_ }
  }

  # combined ranking
  aug$deng_sample <- name
  aug$region      <- REGION
  aug <- aug |>
    dplyr::mutate(
      cm_rank   = dplyr::min_rank(dplyr::desc(score)),
      nb_rank   = dplyr::min_rank(dplyr::desc(replace(nblast, is.na(nblast), -Inf))),
      va_rank   = dplyr::min_rank(dplyr::desc(replace(voxel_attr, is.na(voxel_attr), -Inf))),
      mean_rank = (cm_rank + nb_rank + va_rank) / 3
    ) |>
    dplyr::rename(cm_score = score, cm_n_match = n_match, cm_dx = dx,
                  cm_dy = dy, cm_mirror = mirror)

  out <- dplyr::select(aug, deng_sample, region, banc_id, root_888, preferred_id,
                       super_cluster, cell_function, neurotransmitter_predicted,
                       cm_score, cm_n_match, cm_dx, cm_dy, cm_mirror,
                       nblast, voxel_attr,
                       cm_rank, nb_rank, va_rank, mean_rank) |>
    dplyr::arrange(mean_rank)
  readr::write_csv(out, out_csv)
  cat(sprintf("  wrote %s (%d hits, %d NBLAST, %d voxel_attr)\n",
              out_csv, nrow(out), sum(!is.na(out$nblast)),
              sum(!is.na(out$voxel_attr) & out$voxel_attr > 0)))
  out
}

# --- run all queries -----------------------------------------------
queries <- list.files(QUERY_DIR, pattern = "\\.png$", full.names = TRUE)
cat(sprintf("== %d query MIPs, region '%s' ==\n", length(queries), REGION))
all_hits <- lapply(queries, function(q)
                   tryCatch(process_query(q),
                            error = function(e) { message("FAILED: ", conditionMessage(e)); NULL }))

# --- master CSV ----------------------------------------------------
ok <- !vapply(all_hits, is.null, logical(1))
csvs <- list.files(OUT_DIR, pattern = "_hits\\.csv$", full.names = TRUE)
if (length(csvs)) {
  master <- dplyr::bind_rows(
    lapply(csvs, function(f)
      readr::read_csv(f, show_col_types = FALSE,
                      col_types = readr::cols(banc_id = "c", root_888 = "c",
                                              preferred_id = "c"))))
  master_path <- file.path(OUT_DIR, "master_augmented.csv")
  readr::write_csv(master, master_path)
  cat(sprintf("\nwrote %d rows -> %s\n", nrow(master), master_path))
  print(master |> dplyr::group_by(deng_sample) |>
                  dplyr::summarise(n = dplyr::n(),
                                   n_nblast = sum(!is.na(nblast)),
                                   max_cm   = max(cm_score),
                                   max_nb   = max(nblast, na.rm = TRUE)))
}
