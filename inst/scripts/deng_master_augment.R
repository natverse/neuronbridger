#!/usr/bin/env Rscript
# Augment the per-sample colormip-search CSVs with:
#   - v888 root_888 ID + cell metadata from banc_888_meta.feather
#   - voxel-attribution: each LM voxel votes for its nearest L2-skeleton
#     across the top-K BANC hits. Counts per hit normalised to fraction.
# Re-emits the master CSV.
#
# Inputs:
#   ~/banc_meta/banc_888_meta.feather                           (v888 metadata)
#   ~/banc_meta/banc_banc_space_swc/<root_888>.swc              (L2 skeletons)
#   ~/deng_to_banc/colormip_hits/*_hits.csv                     (per-sample CSVs)
#   ~/deng_to_banc/colormip_search/<name>/elastix/gfp_xform/result.nrrd
#                                                               (warped LM volumes)
#
# Output:
#   ~/deng_to_banc/colormip_hits/master_augmented.csv

suppressMessages({
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Install arrow: install.packages('arrow')")
  library(arrow); library(readr); library(dplyr); library(tibble)
  library(nat); library(reticulate)
  if (!requireNamespace("FNN", quietly = TRUE))
    stop("Install FNN: install.packages('FNN')")
  library(FNN)
})

PYBIN <- "/Users/asbates/Library/r-miniconda-arm64/envs/r-reticulate/bin/python"
reticulate::use_python(PYBIN, required = TRUE)
sitk <- reticulate::import("SimpleITK", convert = FALSE)
np   <- reticulate::import("numpy",     convert = FALSE)

META  <- "~/banc_meta/banc_888_meta.feather"
SWC_DIR <- "~/banc_meta/banc_banc_space_swc"
HITS  <- path.expand("~/deng_to_banc/colormip_hits")
SEARCH_DIR <- path.expand("~/deng_to_banc/colormip_search")

cat("loading v888 metadata...\n")
m888 <- arrow::read_feather(path.expand(META))
m888$root_626 <- as.character(m888$root_626)
m888$root_850 <- as.character(m888$root_850)
m888$root_888 <- as.character(m888$root_888)
cat("metadata rows:", nrow(m888), "\n")

# Augment one sample CSV: add v888 ID + cell metadata + voxel-attribution
augment_sample <- function(csv_path) {
  hits <- readr::read_csv(csv_path, show_col_types = FALSE,
                          col_types = readr::cols(banc_id = "c"))
  if (!nrow(hits)) return(hits)

  # 1. Lookup v888 IDs from root_626 (the colormip library's ID).
  # Per user instruction: prefer root_888 if present, else fall back to root_626.
  lk <- m888 |>
    dplyr::select(root_626, root_888, super_cluster, cell_function,
                  neurotransmitter_predicted) |>
    dplyr::distinct(root_626, .keep_all = TRUE)
  aug <- hits |>
    dplyr::left_join(lk, by = c("banc_id" = "root_626"))
  # canonical id col: root_888 if available, else root_626 (= banc_id)
  aug$preferred_id <- ifelse(is.na(aug$root_888), aug$banc_id, aug$root_888)
  aug$id_source    <- ifelse(is.na(aug$root_888), "root_626", "root_888")

  # 2. Voxel attribution: take each LM fg voxel's BANC nm coord, find its
  # nearest L2-skeleton point across the top-K SWCs, increment that hit's count
  sample_name <- aug$deng_sample[1]
  warp_nrrd <- file.path(SEARCH_DIR, sample_name, "elastix/gfp_xform/result.nrrd")
  if (!file.exists(warp_nrrd)) {
    warning("missing warp NRRD for ", sample_name, "; skipping voxel-attr")
    aug$voxel_attr <- NA_real_
    return(aug)
  }
  img <- sitk$ReadImage(warp_nrrd)
  arr <- np$clip(sitk$GetArrayFromImage(img), 0L, 255L)$astype(np$uint8)
  sp  <- unlist(reticulate::py_to_r(img$GetSpacing()))
  zyx <- np$nonzero(np$greater(arr, np$uint8(20)))
  zs <- reticulate::py_to_r(np$asarray(zyx[[0]]))
  ys <- reticulate::py_to_r(np$asarray(zyx[[1]]))
  xs <- reticulate::py_to_r(np$asarray(zyx[[2]]))
  if (!length(xs)) { aug$voxel_attr <- NA_real_; return(aug) }
  set.seed(7); ix <- sample.int(length(xs), min(length(xs), 25000L))
  pts_um <- cbind(xs[ix]*sp[1], ys[ix]*sp[2], zs[ix]*sp[3])

  # forward-warp via tpsreg (this CSV's region drives which tpsreg)
  region <- aug$region[1]
  if (region == "VNC") {
    suppressMessages(library(bancr))
    library(Morpho)
    tps <- bancr::jrcvnc2018f_to_banc_tpsreg
  } else {
    suppressMessages(library(bancr))
    library(Morpho)
    tps <- bancr::jrc2018f_to_banc_tpsreg
  }
  pts_nm <- Morpho::tps3d(pts_um, refmat = tps$refmat, tarmat = tps$tarmat,
                          lambda = tps$lambda, threads = 4)
  cat(sprintf("  warped %d LM points -> BANC nm\n", nrow(pts_nm)))

  # Read SWCs for top hits — use root_888 if available, else root_626
  swc_paths <- file.path(path.expand(SWC_DIR),
                         paste0(ifelse(is.na(aug$root_888), aug$banc_id, aug$root_888),
                                "_l2.swc"))
  swc_present <- file.exists(swc_paths)
  cat(sprintf("  %d / %d SWCs found locally\n",
              sum(swc_present), length(swc_paths)))
  if (!any(swc_present)) {
    aug$voxel_attr <- NA_real_; return(aug)
  }

  # Build a single point cloud + label vector
  all_pts <- list(); all_lab <- list()
  for (i in which(swc_present)) {
    sk <- tryCatch(nat::read.neuron(swc_paths[i]), error = function(e) NULL)
    if (is.null(sk)) next
    p <- nat::xyzmatrix(sk)   # nm in BANC space (per the bucket convention)
    if (!nrow(p)) next
    all_pts[[length(all_pts)+1]] <- p
    all_lab[[length(all_lab)+1]] <- rep(i, nrow(p))
  }
  if (!length(all_pts)) { aug$voxel_attr <- NA_real_; return(aug) }
  P <- do.call(rbind, all_pts)
  L <- unlist(all_lab)
  cat(sprintf("  built skeleton point cloud: %s pts across %d hits\n",
              format(nrow(P), big.mark=","), length(unique(L))))

  # nearest-neighbour from each LM point
  nn <- FNN::get.knnx(P, pts_nm, k = 1)
  votes <- L[nn$nn.index[,1]]
  # Fraction of total LM points each hit captured
  attr_frac <- tabulate(votes, nbins = nrow(aug)) / length(votes)
  aug$voxel_attr <- attr_frac

  aug
}

# Iterate all CSVs
files <- list.files(HITS, pattern = "_hits\\.csv$", full.names = TRUE)
files <- setdiff(files, c(file.path(HITS, "master.csv"),
                          file.path(HITS, "master_augmented.csv")))
cat("found", length(files), "per-sample CSVs\n")

aug_all <- lapply(files, function(f) {
  cat("\n>", basename(f), "\n")
  augment_sample(f)
})
master <- dplyr::bind_rows(aug_all)

# Combined ranking: average rank across cm_score, nblast (bigger=better),
# voxel_attr (bigger=better). Within each deng_sample.
master <- master |>
  dplyr::group_by(deng_sample) |>
  dplyr::mutate(
    cm_rank   = dplyr::min_rank(dplyr::desc(cm_score)),
    nb_rank   = dplyr::min_rank(dplyr::desc(replace(nblast, is.na(nblast), -Inf))),
    va_rank   = dplyr::min_rank(dplyr::desc(replace(voxel_attr, is.na(voxel_attr), -Inf))),
    mean_rank = (cm_rank + nb_rank + va_rank) / 3
  ) |>
  dplyr::arrange(deng_sample, mean_rank) |>
  dplyr::ungroup()

out <- path.expand("~/deng_to_banc/colormip_hits/master_augmented.csv")
readr::write_csv(master, out)
cat(sprintf("\nwrote %d rows -> %s\n", nrow(master), out))
print(master |> dplyr::group_by(deng_sample) |>
              dplyr::summarise(n = dplyr::n(),
                               n_v888 = sum(!is.na(root_888)),
                               n_nb   = sum(!is.na(nblast)),
                               n_va   = sum(!is.na(voxel_attr) & voxel_attr > 0),
                               max_mean_rank_neuron = banc_id[which.min(mean_rank)]))
