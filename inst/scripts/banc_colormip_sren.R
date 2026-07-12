#!/usr/bin/env Rscript
# End-to-end reproducer for the search_banc_with_lm vignette:
#
#   raw LSM (partial FOV, abdominal-tip-only)
#       -> ch0=GFP, ch2=NC82 NRRDs               (lsm_to_nrrd.py)
#       -> NC82 binary mask                      (lsm_build_nc82_mask.py)
#       -> PCA pre-rotation onto template axes,
#          GFP CoG anchored to (131, 508, 92) um (lsm_pca_rotate.py)
#       -> resample onto JRC2018VNCF grid        (SimpleITK in-script)
#       -> mask GFP by NC82 neuropil             (lsm_neuropil_mask.py)
#       -> bridge to JRC2018VNCU_HR              (jrcvnc2018f_to_jrcvnc2018u_hr_h5)
#       -> colour-MIP query                      (nrrd_to_mip)
#       -> search BANC efferent MIP library      (colormip_search)
#       -> full elastix JRC->BANC transform      (bancr::banc_to_JRC2018F, method="elastix")
#       -> NBLAST fwd+rev mean against full BANC efferent library
#       -> top-6 cMIP + NBLAST panels +
#          side-by-side SREN-vs-top-hit cMIP
#
# Top hits under this pipeline are dominated by BANC cell_type EFFabg40
# (MANC bridge: MNad17), part of the abd_Dsx_OA / abdomen_neurosecretory_cell
# family -- doublesex-expressing sex-dimorphic abdominal neurosecretory
# cells. The hypothesised SREN pair (manc_cell_type EN00B016:
# 720575941680133053 left + 720575941659397456 right) sits in the top
# 10% of the ~1000-cell library across every configuration; visualised
# separately as banc_colormip_sren_en00b016_pair.png.
#
# Registration note: earlier revisions of this script used the bancr
# jrcvnc2018f_to_banc_tpsreg (a 5710-landmark TPS approximation of the
# JRC2018F->BANC alignment) which collapses bilateral signal at the
# abdominal tip. We now use bancr's full elastix chain
# (registrations/vnc_240721: 0_manual_affine + 1_elastix_affine +
# 2_elastix_Bspline_coarse + 3_elastix_Bspline_fine) via
# bancr::banc_to_JRC2018F(..., method="elastix", inverse=TRUE), which
# preserves the bilateral signal correctly. Cell-type annotations come
# from the live SeaTable banc_meta (needs BANCTABLE_TOKEN) rather than
# the frozen banc_888_meta.feather, so both current `cell_type` and
# the `manc_cell_type` bridge are up to date.
#
# Edit SOURCE/OUT_DIR/etc below and run from the repo root.

suppressMessages({
  library(neuronbridger); library(bancr)
  library(arrow); library(dplyr); library(nat); library(nat.ggplot)
  library(ggplot2); library(patchwork); library(scales); library(reticulate)
  library(nat.nblast)
})

# ----- 1. paths -------------------------------------------------------
SOURCE_LSM <- file.path("~/Library/CloudStorage/Dropbox-HMS/Alexander Bates",
                        "neuroanat/KDRC/20250306_C24-53_EN1-5-1060_female_vnc_40x-2.lsm")
SOURCE_LSM <- path.expand(SOURCE_LSM)
WORK    <- "/tmp/banc_colormip_work"
REG_DIR <- file.path(WORK, "sren_register/pca_rot_v2")
LIB     <- file.path(WORK, "library", "mips")
TEMPLATE <- path.expand("~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/templates/JRC2018_VNC_FEMALE_461.nrrd")
OUT_DIR <- "inst/images"

# Target anchor in JRC2018VNCF um -- GFP CoG of the upstream KDRC-
# registered TIF, taken as a known-good rough placement of the SREN
# arborisation in template coords.
TARGET_UM <- "131,508,92"

# LSM channel indices (0-based). For the SREN female LSM: ch0 = GFP,
# ch1 = unused red, ch2 = NC82.
GFP_CHANNEL  <- 0L
NC82_CHANNEL <- 2L

# Hypothesised SREN pair (manc_cell_type EN00B016). Rendered as a
# separate side-by-side figure alongside the actual cMIP-top hit.
SREN_PAIR <- c("720575941680133053",   # left  (BANC cell_type EFFabg09,
                                        #        MANC EN00B016, ACh)
               "720575941659397456")    # right (BANC cell_type EN00B016,
                                        #        MANC EN00B016, octopamine)

META_GS <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather"
MIPS_GS <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_colormips/template_alignment_240721/JRC2018_VNC_UNISEX_461/"

dir.create(WORK,    showWarnings = FALSE, recursive = TRUE)
dir.create(LIB,     showWarnings = FALSE, recursive = TRUE)
dir.create(REG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

py     <- reticulate::conda_python("r-reticulate")
pydir  <- system.file("python", package = "neuronbridger")
py_lsm_to_nrrd        <- file.path(pydir, "lsm_to_nrrd.py")
py_build_mask         <- file.path(pydir, "lsm_build_nc82_mask.py")
py_pca_rotate         <- file.path(pydir, "lsm_pca_rotate.py")
py_neuropil_mask      <- file.path(pydir, "lsm_neuropil_mask.py")
stopifnot(all(file.exists(c(py_lsm_to_nrrd, py_build_mask,
                            py_pca_rotate, py_neuropil_mask))))

# ----- 2. ensure JRC VNC H5 transforms are downloaded -----------------
if (!file.exists(path.expand("~/flybrain-data/JRCVNC2018U_JRCVNC2018F.h5"))) {
  reticulate::py_run_string(
    "import flybrains; flybrains.download_jrc_vnc_transforms()")
}

# ----- 3. LSM -> per-channel NRRDs at native voxdims ------------------
nc82_nrrd <- file.path(REG_DIR, "sren_female_nc82.nrrd")
gfp_nrrd  <- file.path(REG_DIR, "sren_female_gfp.nrrd")
if (!file.exists(nc82_nrrd) || !file.exists(gfp_nrrd)) {
  system2(py, c(shQuote(py_lsm_to_nrrd),
                shQuote(SOURCE_LSM), shQuote(REG_DIR),
                "--prefix", "sren_female",
                "--gfp-channel",  as.integer(GFP_CHANNEL),
                "--nc82-channel", as.integer(NC82_CHANNEL)))
}

# ----- 4. NC82 mask + PCA pre-rotation onto template grid -------------
mask_nrrd <- file.path(REG_DIR, "sren_female_nc82_mask.nrrd")
if (!file.exists(mask_nrrd))
  system2(py, c(shQuote(py_build_mask), shQuote(nc82_nrrd), shQuote(mask_nrrd)))

# Build a posterior-third fixed mask for PCA's default-target path. We
# override the placement with --target-um anyway, but lsm_pca_rotate.py
# still requires a fixed-mask NRRD argument.
fixed_mask <- file.path(REG_DIR, "jrcvnc2018f_tip_mask.nrrd")
if (!file.exists(fixed_mask)) {
  reticulate::py_run_string(sprintf(
    "import SimpleITK as sitk, numpy as np
tpl = sitk.ReadImage('%s')
arr = sitk.GetArrayFromImage(tpl); sx, sy, sz = tpl.GetSpacing()
ny = arr.shape[1]
yidx = ((np.arange(ny) + 0.5) * sy)
m = np.broadcast_to((yidx > 350)[None,:,None], arr.shape).astype(np.uint8).copy()
out = sitk.GetImageFromArray(m); out.CopyInformation(tpl)
sitk.WriteImage(out, '%s', useCompression=True)", TEMPLATE, fixed_mask))
}

rot_dir <- file.path(REG_DIR, "rotated")
if (!file.exists(file.path(rot_dir, "gfp_rot.nrrd"))) {
  dir.create(rot_dir, showWarnings = FALSE, recursive = TRUE)
  system2(py, c(shQuote(py_pca_rotate),
                shQuote(nc82_nrrd), shQuote(gfp_nrrd), shQuote(mask_nrrd),
                shQuote(TEMPLATE),  shQuote(fixed_mask), shQuote(rot_dir),
                "--align-source", "gfp_cog",
                "--target-um", shQuote(TARGET_UM)))
}

# ----- 5. Resample rotated NRRDs onto template grid (uint8 NRRD) ------
gfp_jrcf <- file.path(REG_DIR, "sren_female_gfp_PCAonly_JRC2018VNCF.nrrd")
nc82_jrcf <- file.path(REG_DIR, "sren_female_nc82_PCAonly_JRC2018VNCF.nrrd")
if (!file.exists(gfp_jrcf) || !file.exists(nc82_jrcf)) {
  reticulate::py_run_string(sprintf(
    "import SimpleITK as sitk
tpl = sitk.ReadImage('%s')
for src_in, dst_out in [('%s', '%s'), ('%s', '%s')]:
    s = sitk.ReadImage(src_in)
    ref = sitk.Image(tpl.GetSize(), s.GetPixelID()); ref.CopyInformation(tpl)
    out = sitk.Resample(s, ref, sitk.Transform(), sitk.sitkLinear, 0)
    sitk.WriteImage(out, dst_out, useCompression=True)
    print('wrote', dst_out)",
    TEMPLATE,
    file.path(rot_dir, "gfp_rot.nrrd"), gfp_jrcf,
    file.path(rot_dir, "nc82_rot.nrrd"), nc82_jrcf))
}

# ----- 6. Mask GFP by NC82 neuropil -----------------------------------
gfp_npmask <- file.path(REG_DIR, "sren_female_gfp_PCAonly_neuropilmasked_JRC2018VNCF.nrrd")
if (!file.exists(gfp_npmask))
  system2(py, c(shQuote(py_neuropil_mask),
                shQuote(nc82_jrcf), shQuote(gfp_jrcf), shQuote(gfp_npmask),
                "--mask-out",
                shQuote(file.path(REG_DIR, "neuropil_mask_JRC2018VNCF.nrrd"))))

# Cast uint8 NRRD -> int16 so mmand::morph in nrrd_to_mip accepts it
# (it rejects `raw` storage mode).
recast_int <- function(p_in) {
  v <- nat::read.nrrd(p_in); hdr <- attr(v, "header")
  vox <- abs(diag(hdr[["space directions"]]))
  v_int <- as.integer(v); dim(v_int) <- dim(v)
  nat::im3d(v_int, dims = dim(v), voxdims = vox, origin = c(0, 0, 0))
}
gfp_jrcf_int <- file.path(REG_DIR, "sren_female_in_JRC2018VNCF_int.nrrd")
if (!file.exists(gfp_jrcf_int) || file.mtime(gfp_jrcf_int) < file.mtime(gfp_npmask))
  nat::write.nrrd(recast_int(gfp_npmask), gfp_jrcf_int,
                  dtype = "short", enc = "gzip")

# ----- 7. JRC2018VNCF -> JRC2018VNCU_HR -------------------------------
gfp_hr_raw <- file.path(REG_DIR, "sren_female_in_JRC2018VNCU_HR.nrrd")
gfp_hr     <- file.path(REG_DIR, "sren_female_in_JRC2018VNCU_HR_int.nrrd")
if (!file.exists(gfp_hr_raw) || file.mtime(gfp_hr_raw) < file.mtime(gfp_jrcf_int))
  jrcvnc2018f_to_jrcvnc2018u_hr_h5(input = gfp_jrcf_int, output = gfp_hr_raw,
                                   threads = 8L)
if (!file.exists(gfp_hr) || file.mtime(gfp_hr) < file.mtime(gfp_hr_raw))
  nat::write.nrrd(recast_int(gfp_hr_raw), gfp_hr,
                  dtype = "short", enc = "gzip")

# ----- 8. colour-MIP query --------------------------------------------
query_png <- nrrd_to_mip(
  input        = gfp_hr,
  savefolder   = REG_DIR,
  method       = "direct",
  target_space = "VNC",
  threshold    = 0.80,
  denoise      = "median3d",
  format       = "png",
  save         = TRUE,
  overwrite    = TRUE)
file.copy(query_png,
          file.path(OUT_DIR, "banc_colormip_sren_query.png"),
          overwrite = TRUE)

# ----- 9. Library + cMIP search ---------------------------------------
# Pull annotations from the LIVE SeaTable banc_meta (needs BANCTABLE_TOKEN
# in ~/.Renviron), so we get current cell_type + manc_cell_type. The
# cMIP library filenames + our search results are keyed on root_888
# (frozen v888 IDs), so we join by root_888, not root_id (which drifts
# with further proofreading).
meta_live <- bancr::banctable_query(paste(
  "SELECT root_id, root_888, cell_type, manc_cell_type, cell_class,",
         "cell_sub_class, side, hemilineage, nerve, neuromere,",
         "neurotransmitter_predicted, peripheral_target_type, cell_function,",
         "flow",
  "FROM banc_meta"))
meta <- as.data.frame(meta_live) |>
  filter(flow == "efferent") |>
  select(root_888, side, hemilineage, nerve, neuromere,
         cell_class, cell_sub_class, cell_type, manc_cell_type,
         neurotransmitter_predicted, peripheral_target_type, cell_function)

if (length(list.files(LIB, "\\.png$")) < 700L) {
  all_mips <- system2("gsutil", c("ls", MIPS_GS), stdout = TRUE)
  mip_ids  <- sub(".*/", "",
                  sub("_in_JRC2018_VNC_UNISEX_461\\.png$", "", all_mips))
  to_get   <- all_mips[mip_ids %in% meta$root_888]
  writeLines(to_get, file.path(WORK, "to_get.txt"))
  system(sprintf("xargs -P 16 -I {} gsutil -q cp {} %s/ < %s",
                 LIB, file.path(WORK, "to_get.txt")))
}
cat("library size:", length(list.files(LIB, "\\.png$")), "PNGs\n")

res <- colormip_search(
  query       = query_png,
  library     = LIB,
  threshold   = 100L,
  z_tolerance = 8L, xy_shift = 3L,
  mirror      = FALSE,
  mc.cores    = 8L)
res$root_888 <- sub(".*/(\\d+)_in_JRC2018.*", "\\1", res$path)
saveRDS(res, file.path(REG_DIR, "search_result_npmask.rds"))
cat("\nTop 10 cMIP:\n")
print(head(res[, c("root_888", "score", "n_match", "mirror")], 10))

# ----- 10. SREN dotprops in BANC um (from npmasked GFP) ---------------
# Median-filter the npmasked GFP to suppress isolated noise voxels, but
# keep ALL connected components above the intensity floor -- do NOT
# select the largest-CC only. The previous revision picked the largest
# CC, which dropped the contralateral peak of a genuinely-bilateral
# SREN signal and produced a unilateral query.
ser <- nat::read.nrrd(gfp_jrcf_int)
ser <- mmand::medianFilter(ser, mmand::shapeKernel(c(3,3,3), type = "box"))
storage.mode(ser) <- "integer"
cut <- as.integer(stats::quantile(ser[ser > 0], 0.50))
mask <- ser >= cut
# Drop tiny components (isolated noise) but keep every component with
# >= 50 voxels so both bilateral lobes survive if present.
cc <- mmand::components(mask, mmand::shapeKernel(c(3,3,3), type = "box"))
tab <- tabulate(cc[!is.na(cc)])
mask <- !is.na(cc) & cc %in% which(tab >= 50L)
idx <- which(mask, arr.ind = TRUE)
ser_full <- data.frame(
  X = (idx[,1] - 1L) * 0.461122 + 0.461122/2,
  Y = (idx[,2] - 1L) * 0.461122 + 0.461122/2,
  Z = (idx[,3] - 1L) * 0.7      + 0.35,
  I = as.integer(ser[mask]))
saveRDS(ser_full, file.path(REG_DIR, "sren_pts_jrcvncf_npmask.rds"))
cat("SREN voxels (neuropil-masked, cleaned, CC>=50):", nrow(ser_full), "\n")

# ----- 11. NBLAST fwd+rev mean against the full efferent library ------
# Bridge SREN points JRC2018VNCF -> BANC um using the FULL ELASTIX
# chain (bancr::banc_to_JRC2018F method="elastix" inverse=TRUE), NOT the
# older TPS approximation `jrcvnc2018f_to_banc_tpsreg`. The TPS collapses
# bilateral signal at the abdominal tip; the full elastix chain
# (0_manual_affine + 1_elastix_affine + 2_elastix_Bspline_coarse +
# 3_elastix_Bspline_fine) preserves it.
banc_dps <- if (file.exists(file.path(WORK, "banc_l2dp_efferent_full.rds"))) {
  readRDS(file.path(WORK, "banc_l2dp_efferent_full.rds"))
} else {
  ids_eff <- as.character(meta$root_888)
  cat("Fetching L2 dotprops for", length(ids_eff), "efferents...\n")
  dp <- bancr::banc_read_l2dp(ids_eff, OmitFailures = TRUE)
  saveRDS(dp, file.path(WORK, "banc_l2dp_efferent_full.rds"))
  dp
}
banc_dps <- banc_dps[!sapply(banc_dps, is.null)]
banc_dps <- banc_dps[!duplicated(names(banc_dps))]
cat("BANC L2 dotprops:", length(banc_dps), "\n")

set.seed(11L); N_NB <- 50000L
keep <- sample.int(nrow(ser_full), min(N_NB, nrow(ser_full)),
                   prob = ser_full$I / sum(ser_full$I))
pts_F <- as.matrix(ser_full[keep, c("X","Y","Z")])
pts_BANC_um <- bancr::banc_to_JRC2018F(
  pts_F, region = "vnc", banc.units = "um",
  inverse = TRUE, method = "elastix")
ser_dp <- nat::dotprops(pts_BANC_um, k = 5L)

fwd <- nat.nblast::nblast(ser_dp, nat::as.neuronlist(banc_dps),
                          smat = nat.nblast::smat.fcwb, normalised = TRUE)
rev <- nat.nblast::nblast(nat::as.neuronlist(banc_dps), ser_dp,
                          smat = nat.nblast::smat.fcwb, normalised = TRUE)
nb <- data.frame(root_888   = names(fwd),
                 nblast_fwd = as.numeric(fwd),
                 nblast_rev = as.numeric(rev),
                 nblast_mean = (as.numeric(fwd) + as.numeric(rev)) / 2)
nb_meta <- merge(nb, meta, by = "root_888", all.x = TRUE, sort = FALSE)
nb_meta <- nb_meta[order(-nb_meta$nblast_mean), ]
saveRDS(nb_meta, file.path(REG_DIR, "nblast_scores_fullElastix_fwdrev.rds"))

cat("\nTop 10 (NBLAST fwd+rev mean, full elastix bridge):\n")
print(head(nb_meta[, c("root_888","nblast_mean","cell_type","manc_cell_type",
                       "cell_sub_class","side","neurotransmitter_predicted")], 10))
cat("\nRank of hypothesised SREN pair (EN00B016):\n")
print(nb_meta[nb_meta$root_888 %in% SREN_PAIR,
              c("root_888","nblast_mean","cell_type","manc_cell_type","side",
                "neurotransmitter_predicted")])

# ----- 12. Top-6 hits panels (cMIP + NBLAST) --------------------------
# Same full-elastix bridge as stage 11, on a smaller sample for viz.
set.seed(42L); N_VIZ <- 25000L
keep_viz <- sample.int(nrow(ser_full), min(N_VIZ, nrow(ser_full)),
                       prob = ser_full$I / sum(ser_full$I))
ser_BANC_viz_um <- bancr::banc_to_JRC2018F(
  as.matrix(ser_full[keep_viz, c("X","Y","Z")]),
  region = "vnc", banc.units = "um",
  inverse = TRUE, method = "elastix")
ser_pts <- data.frame(
  X = ser_BANC_viz_um[,1], Y = ser_BANC_viz_um[,2],
  Z = ser_BANC_viz_um[,3], I = ser_full$I[keep_viz])
ser_pal <- scales::col_numeric(
  palette = c("#0a3d2a","#147a4a","#1f8b1f","#7be37b","#d8ff7f"),
  domain  = range(ser_pts$Z))
ser_pts$col <- ser_pal(ser_pts$Z)
hi <- max(ser_pts$I); ser_pts$alpha <- pmax(0.05, (ser_pts$I/hi)^1.5)
ser_pts <- ser_pts[order(-ser_pts$Z), ]

mk_panel <- function(meshes_um, df, score_col, label) function(i) {
  rid <- df$root_888[i]; fmt <- function(x) ifelse(is.na(x) | x == "", "—", x)
  title <- paste(
    sprintf("%s rank %d  —  score %.3f", label, i, df[[score_col]][i]),
    rid,
    sprintf("%s  ->  MANC %s",
            fmt(df$cell_type[i]), fmt(df$manc_cell_type[i])),
    sprintf("%s / %s",
            fmt(df$cell_class[i]), fmt(df$cell_sub_class[i])),
    sprintf("side: %s  NT: %s",
            fmt(df$side[i]), fmt(df$neurotransmitter_predicted[i])),
    sep = "\n")
  base_theme <- theme_void(base_size = 11) +
    theme(plot.title.position = "plot",
          plot.title  = element_text(hjust = 0.5, size = 9.5, lineheight = 1.05),
          plot.margin = margin(4, 8, 4, 8), legend.position = "none")
  mesh_bb <- do.call(rbind, lapply(meshes_um, function(m)
    apply(nat::xyzmatrix(m), 2, range)))
  xlim <- range(c(mesh_bb[,1], ser_pts$X)) + c(-15, 15)
  ylim <- range(c(mesh_bb[,2], ser_pts$Y)) + c(-15, 15)
  zlim <- range(c(mesh_bb[,3], ser_pts$Z)) + c(-10, 10)
  rot_xz <- matrix(c(1,0,0, 0,0,1, 0,1,0), nrow = 3, byrow = TRUE)
  yx <- ggplot() +
    geom_point(data = ser_pts,
               aes(x = X, y = Y, alpha = alpha, colour = col),
               size = 0.45, stroke = 0, show.legend = FALSE) +
    geom_neuron(meshes_um[[rid]], cols = "#aa44cc", alpha = 0.7) +
    scale_alpha_identity() + scale_color_identity() + scale_fill_identity() +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_y_reverse() + labs(title = title) + base_theme
  xz <- ggplot() +
    geom_point(data = ser_pts,
               aes(x = X, y = Z, alpha = alpha, colour = col),
               size = 0.45, stroke = 0, show.legend = FALSE) +
    geom_neuron(meshes_um[[rid]], cols = "#aa44cc", alpha = 0.7,
                rotation_matrix = rot_xz) +
    scale_alpha_identity() + scale_color_identity() + scale_fill_identity() +
    coord_fixed(xlim = xlim, ylim = zlim, expand = FALSE) + base_theme
  patchwork::wrap_plots(yx, xz, ncol = 1, heights = c(3, 1.4))
}

fetch_with_retry <- function(rid, tries = 4L) {
  for (i in seq_len(tries)) {
    m <- try(bancr::banc_read_neuron_meshes(rid), silent = TRUE)
    if (!inherits(m, "try-error") && length(m) == 1L) return(m[[1]])
    Sys.sleep(2L ^ (i - 1L))
  }
  NULL
}
get_meshes <- function(ids) {
  out <- lapply(ids, function(rid) {
    m <- fetch_with_retry(rid); if (is.null(m)) return(NULL)
    nat::xyzmatrix(m) <- nat::xyzmatrix(m) / 1000; m
  })
  names(out) <- ids
  if (any(vapply(out, is.null, logical(1))))
    stop("missing meshes: ",
         paste(ids[vapply(out, is.null, logical(1))], collapse = ", "))
  out
}

# cMIP panel
top_cm <- head(res, 6L) |> left_join(meta, by = "root_888")
meshes_cm <- get_meshes(top_cm$root_888)
p_cm <- patchwork::wrap_plots(
  lapply(seq_len(nrow(top_cm)), mk_panel(meshes_cm, top_cm, "score", "cMIP")),
  ncol = 3) +
  plot_annotation(
    title    = "Top-6 cMIP hits — SREN PCA + neuropil-mask registration",
    subtitle = paste0("BANC EM hit (purple) + SREN LM (hue=depth, alpha=intensity), ",
                      "BANC native space"))
ggsave(file.path(OUT_DIR, "banc_colormip_sren_top_hits.png"),
       p_cm, width = 13, height = 13, dpi = 130, bg = "white")
cat("wrote", file.path(OUT_DIR, "banc_colormip_sren_top_hits.png"), "\n")

# NBLAST fwd+rev mean panel
top_nb <- head(nb_meta, 6L)
meshes_nb <- get_meshes(top_nb$root_888)
p_nb <- patchwork::wrap_plots(
  lapply(seq_len(nrow(top_nb)),
         mk_panel(meshes_nb, top_nb, "nblast_mean", "NBLAST fwd+rev")), ncol = 3) +
  plot_annotation(
    title    = "Top-6 NBLAST fwd+rev hits — full elastix JRC2018F->BANC bridge",
    subtitle = paste0("BANC EM hit (purple) + SREN LM (hue=depth, alpha=intensity), ",
                      "BANC native space"))
ggsave(file.path(OUT_DIR, "banc_colormip_sren_nblast_top.png"),
       p_nb, width = 13, height = 13, dpi = 130, bg = "white")
cat("wrote", file.path(OUT_DIR, "banc_colormip_sren_nblast_top.png"), "\n")

# EN00B016 pair overlay (the hypothesised SREN pair). Panels the two
# candidate cells against the SREN LM cloud so a visual anatomy check
# is possible in the vignette.
pair_df <- meta[meta$root_888 %in% SREN_PAIR, , drop = FALSE]
pair_df$nblast_mean <- nb_meta$nblast_mean[match(pair_df$root_888, nb_meta$root_888)]
pair_df$nblast_rank <- match(pair_df$root_888, nb_meta$root_888)
pair_df$score <- pair_df$nblast_mean
meshes_pair <- get_meshes(pair_df$root_888)
p_pair <- patchwork::wrap_plots(
  lapply(seq_len(nrow(pair_df)), mk_panel(
    meshes_pair, pair_df, "nblast_mean",
    "hypothesised SREN (NBLAST fwd+rev)")),
  ncol = 2) +
  plot_annotation(
    title    = "Hypothesised SREN pair: manc_cell_type EN00B016",
    subtitle = paste0("BANC EM meshes overlaid on SREN LM cloud ",
                      "(full-elastix bridge, BANC um)."))
ggsave(file.path(OUT_DIR, "banc_colormip_sren_en00b016_pair.png"),
       p_pair, width = 11, height = 8, dpi = 130, bg = "white")
cat("wrote", file.path(OUT_DIR, "banc_colormip_sren_en00b016_pair.png"), "\n")

# ----- 13. Side-by-side: SREN query vs top cMIP hit -------------------
# Uses a tight Y-crop on the abdominal tip so both cMIPs read at the
# same scale. Labels the top hit dynamically from live seatable.
TOP_HIT_ROOT <- as.character(res$root_888[1])
top_meta <- meta[meta$root_888 == TOP_HIT_ROOT, , drop = FALSE]
top_hit_label <- if (nrow(top_meta) > 0) {
  sprintf("%s (%s, MANC %s, %s, %s)",
          TOP_HIT_ROOT,
          ifelse(is.na(top_meta$side), "?", top_meta$side),
          ifelse(is.na(top_meta$manc_cell_type), "?", top_meta$manc_cell_type),
          ifelse(is.na(top_meta$cell_type), "?", top_meta$cell_type),
          ifelse(is.na(top_meta$neurotransmitter_predicted), "?",
                 top_meta$neurotransmitter_predicted))
} else TOP_HIT_ROOT
top_hit_png <- file.path(LIB, paste0(TOP_HIT_ROOT,
                                     "_in_JRC2018_VNC_UNISEX_461.png"))
reticulate::py_run_string(sprintf(
"import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
q = np.asarray(Image.open('%s')); h = np.asarray(Image.open('%s'))
nz = (q[...,:3].sum(-1) > 0) | (h[...,:3].sum(-1) > 0)
ys = np.where(nz)[0]; pad = 30
y0 = max(0, ys.min() - pad); y1 = min(q.shape[0], ys.max() + pad)
qc, hc = q[y0:y1], h[y0:y1]
fig, axs = plt.subplots(1, 2, figsize=(10, 7), facecolor='black',
                        gridspec_kw={'wspace': 0.05})
axs[0].imshow(qc); axs[0].set_facecolor('black')
axs[0].set_title('SREN LM query\\nneuropil-masked GFP, in JRC2018VNCU_HR',
                 fontsize=12, color='white', pad=10); axs[0].axis('off')
axs[1].imshow(hc); axs[1].set_facecolor('black')
axs[1].set_title('Top BANC cMIP hit:\\n%s', fontsize=11, color='white',
                 pad=10); axs[1].axis('off')
fig.suptitle('SREN (Chae et al., APDNC4 2026) vs BANC top cMIP hit  —  '
             'cMIPs at the abdominal tip\\nPCA pre-rotation + neuropil mask '
             '+ full-elastix bridge + NBLAST fwd+rev',
             color='white', fontsize=12, y=0.97)
plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.savefig('%s', dpi=140, facecolor='black',
            bbox_inches='tight', pad_inches=0.25)",
  query_png, top_hit_png, top_hit_label,
  file.path(OUT_DIR, "banc_colormip_sren_vs_top_hit.png")))
cat("wrote", file.path(OUT_DIR, "banc_colormip_sren_vs_top_hit.png"), "\n")
cat("done.\n")
