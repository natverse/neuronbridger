#!/usr/bin/env Rscript
# Build the side-by-side colour-MIP demonstration panel for
# vignettes/colormip_direct_vs_fiji.Rmd and the README.
#
# Renders the BANC v888 AstA1 cells (a left/right peptidergic central-brain
# pair, cell_type "AstA1" in `compiled_data/banc_888/banc_888_meta.feather`
# at gs://lee-lab_brain-and-nerve-cord-fly-connectome) through both
# back-ends of nrrd_to_mip() that we can run locally:
#   * method = "direct" (pure R)
#   * method = "python" (BANC's fanc.render_neurons depth_lut via reticulate)
# The FIJI macro is referenced in the vignette as the canonical original
# but cannot be invoked from a non-interactive build.
#
# Provenance:
#   - Mesh: bancr::banc_read_neuron_meshes() -- public BANC segmentation
#     served from gs://zetta_lee_fly_cns_001_kisuk (no CAVE auth required
#     for read). Mirrored at
#     gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/banc_meshes/
#   - Cell-type lookup: gs://lee-lab_brain-and-nerve-cord-fly-connectome/
#       compiled_data/banc_888/banc_888_meta.feather
#       (cell_type == "AstA1" -> two root_888 IDs, one per side)
#   - Sub-sampled vertex point cloud (both sides, 150k pts each, in nm) is
#     cached at inst/extdata/asta1/banc_asta1_points_nm.rds so the script
#     runs without needing to re-fetch the 50MB drc meshes.

suppressMessages({
  library(neuronbridger)
  library(bancr)
  library(nat)
  library(nat.flybrains)
  library(nat.templatebrains)
  library(abind)
})

OUT_DIR <- "inst/images"
CACHE   <- "inst/extdata/asta1/banc_asta1_points_nm.rds"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Load (or regenerate) the AstA1 point cloud cache ---------------
#
# To regenerate from scratch (re-runs the cloudvolume mesh fetch):
#
#   library(bancr); library(arrow)
#   m <- arrow::read_feather(
#     "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather"
#   )
#   ids <- m$root_888[!is.na(m$cell_type) & m$cell_type == "AstA1"]
#   meshes <- bancr::banc_read_neuron_meshes(ids)
#   pts <- lapply(meshes, function(mm) {
#     xyz <- nat::xyzmatrix(mm)
#     xyz <- xyz[sample(nrow(xyz), min(150000L, nrow(xyz))), ]
#     storage.mode(xyz) <- "integer"
#     xyz
#   })
#   names(pts) <- c("right", "left")
#   saveRDS(pts, CACHE)
pts_nm <- readRDS(CACHE)
cat(sprintf("AstA1 cached point cloud: %s sides, %s pts/side\n",
            length(pts_nm),
            paste(sapply(pts_nm, nrow), collapse = "/")))

# ---- 2) Bridge BANC nm -> JRC2018U microns -----------------------------
# Two-hop chain:
#   BANC (nm) -> JRC2018F via bancr's tpsreg (pure R, ships as data)
#   JRC2018F  -> JRC2018U via nat.jrcbrains saalfeldlab H5 transforms
#                (one-time download via download_saalfeldlab_registrations()).
# JRC2018U is the "Unisex" 20x average that the NeuronBridge ColorMIP
# pipeline lives in; the high-resolution variant `JRC2018U_HR` shares
# the same coordinate system but on a smaller (1210 x 566 x 174,
# 0.519 x 0.519 x 1.0 um) grid that NeuronBridge actually serves MIPs at.
suppressMessages(library(nat.jrcbrains))
nat.jrcbrains::register_saalfeldlab_registrations()

all_nm <- do.call(rbind, pts_nm)
cat("Bridging", nrow(all_nm), "points BANC -> JRC2018F -> JRC2018U ...\n")
pts_jrcF <- bancr::banc_to_JRC2018F(all_nm, method = "tpsreg",
                                    banc.units = "nm", region = "brain")
pts_jrc  <- nat.templatebrains::xform_brain(pts_jrcF,
                                            sample = "JRC2018F",
                                            reference = "JRC2018U")
cat("JRC2018U bbox (microns):\n"); print(apply(pts_jrc, 2, range))

# ---- 3) Voxelise into the NeuronBridge JRC2018U_HR grid ----------------
# NeuronBridge serves brain MIPs at JRC2018U_HR: (1210, 566, 174) at
# (0.5189, 0.5189, 1.0) um. We declare that grid as a templatebrain so
# nat::as.im3d puts our points into the same voxel coords as Janelia's
# pre-computed colour-MIPs.
JRC2018U_HR <- nat.templatebrains::templatebrain(
  "JRC2018U_HR",
  dims    = c(1210L, 566L, 174L),
  voxdims = c(0.5189, 0.5189, 1.0),
  units   = "microns"
)
vol <- nat::as.im3d(pts_jrc, JRC2018U_HR)
storage.mode(vol) <- "integer"
cat("Volume dims:", dim(vol), "  occupied voxels:", sum(vol > 0), "\n")

# ---- 4) Run both back-ends ---------------------------------------------
cat("Rendering method = 'direct' ...\n")
mip_r  <- nrrd_to_mip(vol, save = FALSE, method = "direct",
                      target_space = "brain")
cat("Rendering method = 'python' (reticulate -> BANC depth_lut) ...\n")
mip_py <- nrrd_to_mip(vol, save = FALSE, method = "python",
                      target_space = "brain")

cat(sprintf("R   dim: %s\n", paste(dim(mip_r),  collapse = " x ")))
cat(sprintf("py  dim: %s\n", paste(dim(mip_py), collapse = " x ")))
cat(sprintf("Max abs diff = %.4f (= %d / 255 RGB units)\n",
            max(abs(mip_r - mip_py)),
            as.integer(round(max(abs(mip_r - mip_py)) * 255))))
cat(sprintf("Pixels exactly equal: %d / %d (%.4f%%)\n",
            sum(mip_r == mip_py), length(mip_r),
            100 * sum(mip_r == mip_py) / length(mip_r)))

# ---- 5) Single-method demo for the README ------------------------------
png::writePNG(mip_r, file.path(OUT_DIR, "colormip_direct_demo.png"))
cat("wrote: ", file.path(OUT_DIR, "colormip_direct_demo.png"), "\n")

# ---- 6) Side-by-side panel: direct (R) | python (BANC) | abs diff x50 --
diff_amp <- abs(mip_r - mip_py) * 50
diff_amp[diff_amp > 1] <- 1
dim(diff_amp) <- dim(mip_r)
gutter <- array(0, dim = c(8L, dim(mip_r)[2], 3L))
panel  <- abind::abind(mip_r, gutter, mip_py, gutter, diff_amp, along = 1)
panel[panel > 1] <- 1; panel[panel < 0] <- 0
png::writePNG(panel, file.path(OUT_DIR, "colormip_methods_panel.png"))
cat("wrote: ", file.path(OUT_DIR, "colormip_methods_panel.png"), "\n")
cat("Panel dims (y, x, c): ", paste(dim(panel), collapse = " x "), "\n")
