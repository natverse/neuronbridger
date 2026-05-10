#!/usr/bin/env Rscript
# NBLAST a Deng LM volume's top-K hits against the LOCAL L2 SWCs
# (skipping the rate-limited CAVE backend mesh-fetch path).
#
# Args: <warp_nrrd> <hits_csv> <region (brain|VNC)>

suppressMessages({
  devtools::load_all("/Users/asbates/Projects/natverse/neuronbridger", quiet = TRUE)
  library(bancr); library(arrow); library(nat); library(nat.flybrains)
  library(nat.templatebrains); library(nat.jrcbrains); library(nat.nblast)
  nat.jrcbrains::register_saalfeldlab_registrations()
  library(readr); library(tibble); library(reticulate)
})

PYBIN <- "/Users/asbates/Library/r-miniconda-arm64/envs/r-reticulate/bin/python"
reticulate::use_python(PYBIN, required = TRUE)
sitk <- reticulate::import("SimpleITK", convert = FALSE)
np   <- reticulate::import("numpy",     convert = FALSE)

SWC_DIR <- path.expand("~/banc_meta/banc_banc_space_swc")
META    <- "~/banc_meta/banc_888_meta.feather"

args <- commandArgs(trailingOnly = TRUE)
warp_nrrd <- args[1]
csv_path  <- args[2]
region    <- args[3]
stopifnot(file.exists(warp_nrrd), file.exists(csv_path))

hits <- readr::read_csv(csv_path, show_col_types = FALSE,
                        col_types = readr::cols(banc_id = "c"))

# Lookup root_888 from metadata
m888 <- arrow::read_feather(path.expand(META))
m888$root_626 <- as.character(m888$root_626); m888$root_888 <- as.character(m888$root_888)
lk <- dplyr::distinct(dplyr::select(m888, root_626, root_888), root_626, .keep_all = TRUE)
hits <- dplyr::left_join(hits, lk, by = c("banc_id" = "root_626"))
hits$preferred_id <- ifelse(is.na(hits$root_888), hits$banc_id, hits$root_888)

# Build LM dotprops
img <- sitk$ReadImage(warp_nrrd)
arr <- np$clip(sitk$GetArrayFromImage(img), 0L, 255L)$astype(np$uint8)
sp  <- unlist(reticulate::py_to_r(img$GetSpacing()))
zyx <- np$nonzero(np$greater(arr, np$uint8(20)))
zs <- reticulate::py_to_r(np$asarray(zyx[[0]]))
ys <- reticulate::py_to_r(np$asarray(zyx[[1]]))
xs <- reticulate::py_to_r(np$asarray(zyx[[2]]))
n <- length(xs)
if (!n) { cat("LM has no fg voxels\n"); quit() }
set.seed(1); idx <- if (n > 50000L) sample.int(n, 50000L) else seq_len(n)
lm_pts <- cbind(xs[idx]*sp[1], ys[idx]*sp[2], zs[idx]*sp[3])
lm_dp  <- nat::dotprops(lm_pts, k = 5L)
cat("LM dotprops:", length(lm_dp$points), "points\n")

# Per-hit: load SWC, warp BANC nm -> JRC2018F um (or JRCVNC2018F)
# -> JRC2018U / JRCVNC2018U (matches LM frame), dotprops, NBLAST.
target_tb  <- if (region == "brain") "JRC2018U"  else "JRCVNC2018U"
source_tb  <- if (region == "brain") "JRC2018F"  else "JRCVNC2018F"
region_tps <- if (region == "brain") "brain" else "vnc"

scores <- vapply(seq_len(nrow(hits)), function(i) {
  swc_path <- file.path(SWC_DIR, paste0(hits$preferred_id[i], "_l2.swc"))
  if (!file.exists(swc_path)) return(NA_real_)
  tryCatch({
    sk <- nat::read.neuron(swc_path)
    pts_nm <- nat::xyzmatrix(sk)        # BANC nm
    if (nrow(pts_nm) < 6L) return(NA_real_)
    # BANC nm -> JRC2018F/JRCVNC2018F um via bancr tps
    pts_jrc <- bancr::banc_to_JRC2018F(pts_nm,
                                        method = "tpsreg",
                                        banc.units = "nm",
                                        region = region_tps)
    # ... -> target template (JRC2018U / JRCVNC2018U) using saalfeldlab
    pts_tg <- nat.templatebrains::xform_brain(pts_jrc,
                                              sample = source_tb,
                                              reference = target_tb)
    if (nrow(pts_tg) < 6L) return(NA_real_)
    dp <- nat::dotprops(pts_tg, k = 5L)
    as.numeric(nat.nblast::nblast(dp, target = lm_dp, normalised = TRUE))
  }, error = function(e) NA_real_)
}, numeric(1))

# Only fill in NAs (don't regress prior scores)
prior <- if ("nblast" %in% names(hits)) hits$nblast else rep(NA_real_, nrow(hits))
new <- prior
fill_idx <- which(is.na(prior) & !is.na(scores))
new[fill_idx] <- scores[fill_idx]
hits$nblast <- new

# Drop helper cols before writing back
hits <- dplyr::select(hits, -dplyr::any_of(c("root_888", "preferred_id")))
readr::write_csv(hits, csv_path)
cat(sprintf("updated %s : %d / %d have NBLAST (added %d this run)\n",
            csv_path, sum(!is.na(new)), length(new), length(fill_idx)))
