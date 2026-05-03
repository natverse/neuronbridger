#!/usr/bin/env Rscript
# Build the side-by-side colour-MIP demonstration panel for
# vignettes/colormip_direct_vs_fiji.Rmd and the README.
#
# Renders a synthetic-but-recognisable JRC2018U_HR-shaped 3D volume (a
# branching, depth-varying "pseudo-neuron") through both back-ends of
# nrrd_to_mip() that we can run locally:
#   * method = "direct" (pure R)
#   * method = "python" (BANC's fanc.render_neurons depth_lut via reticulate)
# The FIJI macro is referenced in the vignette as the canonical original
# but cannot be invoked from a non-interactive build, so the panel is
# captioned as "the same algorithm, two implementations".

suppressMessages({
  library(neuronbridger)
  library(abind)
})

OUT_DIR <- "inst/images"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Synthetic JRC2018U_HR volume with a recognisable depth-varying path
nx <- 1210L; ny <- 566L; nz <- 174L
vol <- array(0L, dim = c(nx, ny, nz))

# Anterior dendritic tuft (small) -> long descending arbour -> posterior tuft.
draw_blob <- function(vol, cx, cy, cz, r) {
  xr <- max(1L, cx - r):min(dim(vol)[1], cx + r)
  yr <- max(1L, cy - r):min(dim(vol)[2], cy + r)
  zr <- max(1L, cz - 1L):min(dim(vol)[3], cz + 1L)
  for (x in xr) for (y in yr) {
    if ((x - cx)^2 + (y - cy)^2 <= r^2)
      vol[x, y, zr] <<- 255L
  }
}

n_pts <- 1200L
ts <- seq(0, 1, length.out = n_pts)
for (i in seq_along(ts)) {
  t <- ts[i]
  cx <- as.integer(180 + (nx - 360) * t)
  cy <- as.integer(ny / 2 + 110 * sin(t * 4 * pi) * (1 - t * 0.3))
  cz <- as.integer(8 + (nz - 16) * t)
  if (cx >= 1 && cx <= nx && cy >= 1 && cy <= ny && cz >= 1 && cz <= nz)
    vol[cx, cy, cz] <- 255L
  if (cx + 1 <= nx) vol[cx + 1, cy, cz] <- 255L
  if (cy + 1 <= ny) vol[cx, cy + 1, cz] <- 255L
}
# Tufts at both ends
draw_blob(vol, 200,  ny %/% 2,  10,  18)
draw_blob(vol, nx - 200, ny %/% 2, nz - 10, 22)
# Two side branches
for (i in seq_along(ts)) {
  t <- ts[i]
  if (t > 0.30 && t < 0.34) {
    cx <- as.integer(180 + (nx - 360) * t)
    cy <- as.integer(ny / 2 + 200 * (t - 0.30) / 0.04)
    cz <- as.integer(8 + (nz - 16) * t)
    if (cx <= nx && cy <= ny && cz <= nz) vol[cx, cy, cz] <- 255L
  }
  if (t > 0.62 && t < 0.66) {
    cx <- as.integer(180 + (nx - 360) * t)
    cy <- as.integer(ny / 2 - 200 * (t - 0.62) / 0.04)
    cz <- as.integer(8 + (nz - 16) * t)
    if (cx <= nx && cy >= 1 && cz <= nz) vol[cx, cy, cz] <- 255L
  }
}

# ---- 2) Run both back-ends
cat("Rendering method = 'direct' ...\n")
mip_r <- nrrd_to_mip(vol, save = FALSE, method = "direct",
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

# ---- 3) Single-method demo for the README
png::writePNG(mip_r, file.path(OUT_DIR, "colormip_direct_demo.png"))
cat("wrote: ", file.path(OUT_DIR, "colormip_direct_demo.png"), "\n")

# ---- 4) Side-by-side panel: direct (R) | python (BANC) | abs diff x50
diff_amp <- abs(mip_r - mip_py) * 50           # amplify for visibility
diff_amp[diff_amp > 1] <- 1
dim(diff_amp) <- dim(mip_r)                     # pmin / arithmetic strip dim

# Stack the three RGB images vertically with thin black gutter rows
# (no magick headers — environments without system fonts can't annotate).
gutter <- array(0, dim = c(8L, dim(mip_r)[2], 3L))
panel  <- abind::abind(mip_r, gutter, mip_py, gutter, diff_amp, along = 1)
panel[panel > 1] <- 1; panel[panel < 0] <- 0
png::writePNG(panel, file.path(OUT_DIR, "colormip_methods_panel.png"))
cat("wrote: ", file.path(OUT_DIR, "colormip_methods_panel.png"), "\n")
cat("Panel dims (y, x, c): ", paste(dim(panel), collapse = " x "), "\n")
