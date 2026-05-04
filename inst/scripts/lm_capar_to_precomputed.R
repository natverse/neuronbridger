#!/usr/bin/env Rscript
# Convert a Kondo-2020 receptor-expression NRRD (IS2 template space) to a
# Neuroglancer "precomputed" image layer — the same on-disk format BANC's
# `JRC2018F_brain.ng` and `JRC2018F_aligned240721_to_BANC.ng` atlas layers
# use (browseable via spelunker).
#
# This script demonstrates the precomputed conversion only. To overlay the
# layer on BANC EM in spelunker the volume should first be brought into
# JRC2018F space (or directly into BANC voxel space) — see
# vignettes/lm_layer_neuroglancer.Rmd for the full pipeline.
#
# Source NRRD lives in a private Dropbox path; this script also subsamples
# (4x xy, 2x z) and 8-bit-truncates so the demo output is small (~7 MB).

suppressMessages({
  library(neuronbridger)
  library(nat)
  library(reticulate)
})

NRRD_PATH <- "/Users/asbates/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/kondo_et_al_2020/nrrd/IS2_CapaR_no1_02_warp_m0g40c4e1e-1x16r3.nrrd"
OUT_DIR   <- "/tmp/CapaR_IS2_pc"

stopifnot(file.exists(NRRD_PATH))
v   <- nat::read.nrrd(NRRD_PATH)
hdr <- attr(v, "header")
src_res_um <- diag(hdr[["space directions"]])

cat("source dim    :", dim(v), "\n")
cat("source res um :", round(src_res_um, 4), "\n")

# Down-sample 4x xy / 2x z and squash 16-bit -> 8-bit.
v8 <- v[ seq(1L, dim(v)[1], by = 4L),
         seq(1L, dim(v)[2], by = 4L),
         seq(1L, dim(v)[3], by = 2L) ]
v8 <- as.integer(pmin(pmax(v8, 0L), 4095L) / 16L)
dim(v8) <- c(192L, 192L, 87L)

cat("output dim    :", dim(v8), "\n")
cat("output res nm :", round(src_res_um * c(4, 4, 2) * 1000), "\n")

nrrd_to_precomputed(
  input      = v8,
  output     = OUT_DIR,
  resolution = round(src_res_um * c(4, 4, 2) * 1000),
  data_type  = "uint8",
  encoding   = "raw",
  chunk_size = c(64L, 64L, 32L),
  overwrite  = TRUE
)

# Total size:
files <- list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)
size_mb <- sum(file.size(files)) / 1024^2
cat(sprintf("\nWrote %d files, total %.1f MB.\n", length(files), size_mb))
cat("Upload (lee-lab maintainers only — pick your own bucket otherwise):\n")
cat(sprintf(
  "  gsutil -m cp -r %s gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/kondo_et_al_2020/CapaR_no1_02_IS2/\n",
  OUT_DIR))
