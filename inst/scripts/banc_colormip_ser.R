#!/usr/bin/env Rscript
# End-to-end reproducer for the banc_colormip_search vignette:
#   SER (JRCVNC2018F) -> JRC2018VNCU_HR -> color-MIP -> search vs
#   BANC efferent library -> top-hit nat.ggplot panel.
#
# Run from the repo root after editing SOURCE / WORK below.

suppressMessages({
  library(neuronbridger); library(bancr)
  library(arrow); library(dplyr); library(nat); library(nat.ggplot)
  library(ggplot2); library(patchwork); library(reticulate)
})

# ----- 1. paths --------------------------------------------------------
SOURCE  <- "/path/to/ser_merged registration_female.tif"
WORK    <- "/tmp/banc_colormip_work"
LIB     <- file.path(WORK, "library", "mips")
META_GS <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather"
MIPS_GS <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_colormips/template_alignment_240721/JRC2018_VNC_UNISEX_461/"
OUT_DIR <- "inst/images"
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
dir.create(LIB,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----- 2. ensure JRC VNC H5 transforms are downloaded ------------------
if (!file.exists(path.expand("~/flybrain-data/JRCVNC2018U_JRCVNC2018F.h5"))) {
  reticulate::py_run_string(
    "import flybrains; flybrains.download_jrc_vnc_transforms()")
}

# ----- 3. split TIF channels at JRCVNC2018F native voxdims -------------
nrrds <- tif_to_nrrd_channels(
  input      = SOURCE,
  output_dir = WORK,
  voxdims    = c(0.4, 0.4, 0.4),
  prefix     = "ser_female",
  channels   = c(0L, 1L),
  labels     = c("nc82", "ser"))

# ----- 4. bridge SER channel into JRC2018VNCU_HR -----------------------
ser_in_HR <- jrcvnc2018f_to_jrcvnc2018u_hr_h5(
  input  = nrrds["ser"],
  output = file.path(WORK, "ser_female_in_JRC2018VNCU_HR.nrrd"),
  threads = 8L)

# ----- 5. query color-MIP ----------------------------------------------
query_png <- nrrd_to_mip(
  input        = ser_in_HR,
  savefolder   = WORK,
  method       = "direct",
  target_space = "VNC",
  threshold    = 0.80,
  denoise      = "median3d",
  format       = "png",
  save         = TRUE,
  overwrite    = TRUE)

file.copy(query_png,
          file.path(OUT_DIR, "banc_colormip_ser_query.png"),
          overwrite = TRUE)

# ----- 6. download metadata + efferent MIP library ---------------------
META_TMP <- file.path(WORK, "banc_888_meta.feather")
if (!file.exists(META_TMP))
  system2("gsutil", c("cp", META_GS, META_TMP))

meta <- arrow::read_feather(META_TMP) |>
  filter(flow == "efferent") |>
  select(root_888, side, hemilineage, nerve, neuromere,
         cell_class, cell_type, neurotransmitter_predicted,
         peripheral_target_type, cell_function)

all_mips <- system2("gsutil", c("ls", MIPS_GS), stdout = TRUE)
mip_ids  <- sub(".*/", "",
                sub("_in_JRC2018_VNC_UNISEX_461\\.png$", "", all_mips))
to_get   <- all_mips[mip_ids %in% meta$root_888]
writeLines(to_get, file.path(WORK, "to_get.txt"))
system(sprintf("xargs -P 16 -I {} gsutil -q cp {} %s/ < %s",
               LIB, file.path(WORK, "to_get.txt")))
cat("library size:", length(list.files(LIB, "\\.png$")), "PNGs\n")

# ----- 7. search --------------------------------------------------------
res <- colormip_search(
  query       = query_png,
  library     = LIB,
  threshold   = 100L,
  z_tolerance = 2L,
  xy_shift    = 2L,
  mirror      = TRUE,
  mc.cores    = 8L)
res$root_888 <- sub(".*/(\\d+)_in_JRC2018.*", "\\1", res$path)
saveRDS(res, file.path(WORK, "search_result.rds"))

cat("\nTop 6:\n")
print(head(res[, c("root_888","score","n_match","mirror")], 6))

# ----- 8. nat.ggplot panel: hit (purple) + SER LM (green) in JRC2018VNCF
top    <- head(res, 6) |> left_join(meta, by = "root_888")

# Render in BANC native um (no TPS extrapolation at the posterior tip).
banc_meshes <- bancr::banc_read_neuron_meshes(top$root_888)
meshes_um <- lapply(banc_meshes, function(m) {
  nat::xyzmatrix(m) <- nat::xyzmatrix(m) / 1000
  m
})
names(meshes_um) <- top$root_888

# Build the SER point cloud: median3d denoise -> q=0.50 threshold ->
# largest connected component -> intensity-weighted decimation to 25k
# points -> alpha = (I/max)^2 with a 0.04 floor.
ser <- nat::read.nrrd(nrrds[["ser"]])
ser <- mmand::medianFilter(ser, mmand::shapeKernel(c(3,3,3), type="box"))
storage.mode(ser) <- "integer"

cut  <- as.integer(stats::quantile(ser[ser > 0], 0.50))
mask <- ser >= cut
cc   <- mmand::components(mask, mmand::shapeKernel(c(3,3,3), type="box"))
mask <- !is.na(cc) & cc == which.max(tabulate(cc[!is.na(cc)]))
idx  <- which(mask, arr.ind = TRUE)

ser_full <- data.frame(
  X = (idx[, 1] - 1L) * 0.4 + 0.2,
  Y = (idx[, 2] - 1L) * 0.4 + 0.2,
  Z = (idx[, 3] - 1L) * 0.4 + 0.2,
  I = as.integer(ser[mask]))

set.seed(42L)
N_VIZ <- 25000L
keep <- sample.int(nrow(ser_full), N_VIZ,
                   prob = ser_full$I / sum(ser_full$I))
# JRCVNC2018F um -> BANC um for the panel.
ser_BANC <- Morpho::tps3d(
  as.matrix(ser_full[keep, c("X","Y","Z")]),
  bancr::jrcvnc2018f_to_banc_tpsreg$refmat,
  bancr::jrcvnc2018f_to_banc_tpsreg$tarmat,
  lambda = bancr::jrcvnc2018f_to_banc_tpsreg$lambda)
ser_pts <- data.frame(
  X = ser_BANC[, 1] / 1000, Y = ser_BANC[, 2] / 1000,
  Z = ser_BANC[, 3] / 1000, I = ser_full$I[keep])
hi <- max(ser_pts$I)
ser_pts$alpha <- pmax(0.04, (ser_pts$I / hi) ^ 2)

mesh_bb <- do.call(rbind, lapply(meshes_um, function(m) {
  apply(nat::xyzmatrix(m), 2, range)
}))
xlim <- range(c(mesh_bb[, 1], ser_pts$X)) + c(-15, 15)
ylim <- range(c(mesh_bb[, 2], ser_pts$Y)) + c(-15, 15)

mk_panel <- function(i) {
  rid <- top$root_888[i]
  fmt <- function(x) ifelse(is.na(x), "—", x)
  title <- paste(
    sprintf("rank %d  —  score %.3f", i, top$score[i]),
    rid,
    fmt(top$cell_class[i]),
    fmt(top$cell_type[i]),
    sprintf("NT: %s", fmt(top$neurotransmitter_predicted[i])),
    sep = "\n")
  ggplot() +
    geom_point(data = ser_pts,
               aes(x = X, y = Y, alpha = alpha),
               colour = "#1f8b1f", size = 0.4, show.legend = FALSE) +
    geom_neuron(meshes_um[[rid]], cols = "#aa44cc", alpha = 0.75) +
    scale_alpha_identity() +
    scale_color_identity() + scale_fill_identity() +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_y_reverse() +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 9.5,
                                    lineheight = 1.05),
          legend.position = "none")
}

p <- wrap_plots(lapply(seq_len(nrow(top)), mk_panel), ncol = 3) +
  plot_annotation(
    title    = "Top-6 BANC efferent hits for the SER query (BANC native um)",
    subtitle = paste0("BANC EM hit (purple) + SER LM source (green, ",
                      "alpha = intensity), in BANC native space"))
ggsave(file.path(OUT_DIR, "banc_colormip_ser_top_hits.png"),
       p, width = 13, height = 10, dpi = 130, bg = "white")

cat("colormip panel done.\n")

# ----- 9. Bridge SER NRRD into BANC voxel space (precomputed publish) ---
#
# Uses bancr's vnc_240721 transformix params (BANC_to_template.txt) to
# resample JRCVNC2018F -> BANC voxel grid (2622 x 2950 x 789 @ 0.4 um),
# then casts to uint8 and writes a Neuroglancer precomputed layer to GCS.
#
# Note: cloudvolume reads $GOOGLE_APPLICATION_CREDENTIALS via google.auth.
# If that's set to a path that doesn't exist on this machine (eg. an
# `~/.Renviron` line copied from another box), override it explicitly:
# Sys.setenv(GOOGLE_APPLICATION_CREDENTIALS =
#   "~/.config/gcloud/application_default_credentials.json")

ser_BANC_nrrd <- file.path(WORK, "ser_female_in_BANC_VNC.nrrd")
if (!file.exists(ser_BANC_nrrd)) {
  message("Run transformix with bancr::vnc_240721 (BANC_to_template.txt) ",
          "on ", nrrds[["ser"]], " to produce ", ser_BANC_nrrd,
          " (uint8). See ?bancr::jrc2018f_to_banc_elastix for an R wrapper.")
} else {
  GCS_LM <- paste0("gs://lee-lab_brain-and-nerve-cord-fly-connectome/",
                   "light_level/kdrc/SER_female_aligned240721_to_BANC.ng")
  nrrd_to_precomputed(input = ser_BANC_nrrd, output = GCS_LM,
                      resolution = c(400L, 400L, 400L),
                      data_type = "uint8",
                      encoding  = "raw",
                      chunk_size = c(128L, 128L, 32L),
                      compress = TRUE, overwrite = FALSE)
}

# ----- 10. Per-match Spelunker URLs ------------------------------------
LM_URL <- paste0("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-",
                 "connectome/light_level/kdrc/",
                 "SER_female_aligned240721_to_BANC.ng")
BASE   <- paste0("https://spelunker.cave-explorer.org/#!middleauth+",
                 "https://global.daf-apis.com/nglstate/api/v1/",
                 "6450802162925568")
top$ngl_url <- vapply(top$root_888, function(rid) {
  bancr::banc_lm_scene(lm_url = LM_URL, layer_name = "KDRC SER (female)",
                       range = c(1, 60), opacity = 0.55,
                       blend = "additive",
                       ids = as.character(rid), url = BASE,
                       shorten = TRUE)
}, character(1))
write.csv(top[, c("root_888","score","cell_type",
                  "neurotransmitter_predicted","ngl_url")],
          file.path(WORK, "ngl_links_cmip.csv"), row.names = FALSE)

# ----- 11. NBLAST cross-check on top-100 colorMIP candidates ------------
# Intensity-weighted: sample 50k voxels from ser_full with sampling
# probability proportional to voxel intensity, so bright pixels get
# more representation in the dotprops vector field.
candidates <- head(res$root_888, 100L)
banc_dps   <- bancr::banc_read_l2dp(candidates)

set.seed(11L)
N_NBLAST  <- 50000L
keep_nb   <- sample.int(nrow(ser_full), N_NBLAST,
                        prob = ser_full$I / sum(ser_full$I))
ser_pts_F <- as.matrix(ser_full[keep_nb, c("X","Y","Z")])
ser_BANC_nm <- Morpho::tps3d(
  ser_pts_F,
  bancr::jrcvnc2018f_to_banc_tpsreg$refmat,
  bancr::jrcvnc2018f_to_banc_tpsreg$tarmat,
  lambda = bancr::jrcvnc2018f_to_banc_tpsreg$lambda)
ser_dp_BANC <- nat::dotprops(ser_BANC_nm / 1000, k = 5L)

fwd <- nat.nblast::nblast(ser_dp_BANC, nat::as.neuronlist(banc_dps),
                          smat = nat.nblast::smat.fcwb, normalised = TRUE)
rev <- nat.nblast::nblast(nat::as.neuronlist(banc_dps), ser_dp_BANC,
                          smat = nat.nblast::smat.fcwb, normalised = TRUE)
nb <- data.frame(root_888    = names(banc_dps),
                 nblast_mean = as.numeric((fwd + rev) / 2))
nb <- nb[order(-nb$nblast_mean), ]
saveRDS(nb, file.path(WORK, "nblast_scores.rds"))

cat("done.\n")
