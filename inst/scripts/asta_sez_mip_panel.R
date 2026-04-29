#!/usr/bin/env Rscript
# Side-by-side colour-depth MIP comparison:
#   SS32423 NeuronBridge brain MIP  |  CB0602 (FlyWire) depth-coloured  |  CB0239 depth-coloured
# Candidate panels use nat.ggplot::geom_neuron with the same
# blue→cyan→green→yellow→orange→red→magenta depth ramp NeuronBridge uses,
# so the three panels read consistently.

suppressMessages({
  library(nat); library(nat.ggplot); library(nat.flybrains); library(elmr)
  library(ggplot2); library(grid); library(magick); library(patchwork)
})

OUT_DIR <- "inst/images"
CACHE   <- "inst/extdata/asta_sez/cache"

CB0602  <- "720575940632295751"
CB0239  <- "720575940634548194"

sk <- readRDS(file.path(CACHE, "l2skel_top2.rds"))
names(sk) <- c("CB0602", "CB0239")[match(names(sk), c(CB0602, CB0239))]

brain  <- elmr::FAFB14.surf
np_all <- elmr::FAFB14NP.surf
SEZ    <- intersect(c("GNG","SAD","FLA_R","FLA_L","AMMC_R","AMMC_L"),
                    np_all$RegionList)
sez_surf <- subset(np_all, SEZ)

# NeuronBridge / FlyLight 'PsychedelicRainBow' depth ramp.
DEPTH_LUT <- c("#0033FF","#00CCFF","#00FF66","#CCFF00","#FFCC00",
               "#FF3300","#FF00CC")

# Use brain outline at low alpha to anchor the viewer (matches the
# whole-brain frame of the SS32423 NB MIP); cells get the depth ramp.
mk_depth_panel <- function(neuron, title) {
  ggplot() +
    geom_neuron(brain, cols = "grey25", alpha = 0.20) +
    geom_neuron(neuron, cols = DEPTH_LUT, alpha = 1, lwd = 0.5) +
    coord_fixed() +
    scale_y_reverse() +
    guides(colour = "none", fill = "none") +
    labs(title = title) +
    theme_void(base_size = 14) +
    theme(plot.background     = element_rect(fill = "black", colour = NA),
          panel.background    = element_rect(fill = "black", colour = NA),
          plot.margin         = margin(6, 6, 6, 6),
          plot.title          = element_text(hjust = 0.5, colour = "white",
                                             margin = margin(b = 4)),
          legend.position     = "none")
}

p_cb0602 <- mk_depth_panel(sk[["CB0602"]], "")
p_cb0239 <- mk_depth_panel(sk[["CB0239"]], "")

# Render each candidate at the same aspect ratio as the SS32423 MIP
# (~2.4:1) so they tile cleanly into a 3-row stack.
ggsave(file.path(OUT_DIR, "_tmp_cb0602.png"), p_cb0602,
       width = 14, height = 5.8, dpi = 200, bg = "black")
ggsave(file.path(OUT_DIR, "_tmp_cb0239.png"), p_cb0239,
       width = 14, height = 5.8, dpi = 200, bg = "black")

# 3 rows, each at full target width. Reads top-to-bottom: SS32423 then
# CB0602 then CB0239 — directly comparable.
ss_mip   <- magick::image_read(file.path(OUT_DIR, "asta_sez_SS32423_brain_mip.png"))
cb0602_p <- magick::image_read(file.path(OUT_DIR, "_tmp_cb0602.png"))
cb0239_p <- magick::image_read(file.path(OUT_DIR, "_tmp_cb0239.png"))

target_w <- 2400
ss_mip   <- magick::image_scale(ss_mip,   sprintf("%dx", target_w))
cb0602_p <- magick::image_scale(cb0602_p, sprintf("%dx", target_w))
cb0239_p <- magick::image_scale(cb0239_p, sprintf("%dx", target_w))

mk_hdr <- function(text, h = 70, size = 34) {
  magick::image_blank(target_w, h, color = "black") |>
    magick::image_annotate(text, color = "white", size = size,
                           gravity = "center")
}
panel <- magick::image_append(c(
  mk_hdr("SS32423  —  NeuronBridge brain MIP (v3_9_0, FlyLight Split-GAL4)"),
  ss_mip,
  mk_hdr("CB0602  —  SS32423 ∩ R65D05 (FAFB-v783 L2 skeleton, depth-encoded)", size = 30),
  cb0602_p,
  mk_hdr("CB0239  —  SS32423 ∩ VT019900 (FAFB-v783 L2 skeleton, depth-encoded)", size = 30),
  cb0239_p
), stack = TRUE)
magick::image_write(panel,
                    file.path(OUT_DIR, "asta_sez_mip_panel.png"),
                    format = "png")

# Cleanup tmp files.
file.remove(file.path(OUT_DIR, c("_tmp_cb0602.png", "_tmp_cb0239.png")))
cat("wrote asta_sez_mip_panel.png\n")
