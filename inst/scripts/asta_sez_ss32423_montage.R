#!/usr/bin/env Rscript
# Download ALL 16 SS32423 brain MIPs in NB v3_9_0 and tile them into a
# montage so the reader can see the full set of stainings the pipeline
# is choosing from. NB ships a representative MIP per slide, and the
# 16 brain MIPs of SS32423 are different specimens — useful to inspect
# expression variability across flies.

suppressMessages({
  library(devtools); load_all(".")
  library(magick)
})
NB_VERSION <- "v3_9_0"
OUT_DIR <- "inst/images"
TMP_DIR <- file.path("inst", "extdata", "asta_sez", "cache", "ss_mips")
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

info <- neuronbridge_info("SS32423", dataset = "by_line",
                          version = NB_VERSION)
brain <- info[info$anatomicalArea == "Brain", ]
cat("brain MIP count:", nrow(brain), "\n")

# Download each brain MIP from S3.
files <- character(nrow(brain))
for (i in seq_len(nrow(brain))) {
  url <- paste0("https://s3.amazonaws.com/janelia-flylight-color-depth/",
                brain$files.CDM[i])
  out <- file.path(TMP_DIR, sprintf("SS32423_brain_%02d.png", i))
  if (!file.exists(out)) {
    try(download.file(url, out, mode = "wb", quiet = TRUE), silent = TRUE)
  }
  files[i] <- out
}
files <- files[file.exists(files)]
cat("downloaded:", length(files), "\n")

# Read, scale to a common width, lay out 4 cols × 4 rows.
imgs <- magick::image_read(files)
imgs <- magick::image_scale(imgs, "600x")
montage <- magick::image_montage(imgs, tile = "4x4",
                                 geometry = "600x300+8+8",
                                 bg = "black")
# Footer label.
w <- magick::image_info(montage)$width
hdr <- magick::image_blank(w, 80, color = "black") |>
  magick::image_annotate(
    sprintf("SS32423 — all %d brain MIPs in NeuronBridge v3_9_0", length(files)),
    color = "white", size = 36, gravity = "center")
final <- magick::image_append(c(hdr, montage), stack = TRUE)
magick::image_write(final, file.path(OUT_DIR, "asta_sez_SS32423_montage.png"),
                    format = "png")
cat("wrote asta_sez_SS32423_montage.png  (", length(files), "MIPs)\n")
