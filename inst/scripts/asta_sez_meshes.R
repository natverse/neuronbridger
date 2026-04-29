#!/usr/bin/env Rscript
# Re-render the AstA candidate panels from FlyWire MESHES (via
# fafbseg::read_cloudvolume_meshes) rather than L2 skeletons. Meshes
# carry the volumetric arbour and so produce a fuller 2D MIP-like
# render under nat.ggplot::geom_neuron().
#
# Cost: meshes are ~20 MB each; ~10 s/neuron to fetch the first time.
# Cached to inst/extdata/asta_sez/cache/meshes_*.rds.

suppressMessages({
  library(devtools); load_all(".")
  library(nat); library(nat.ggplot); library(elmr); library(fafbseg)
  library(ggplot2); library(patchwork); library(magick)
})

OUT_DIR <- "inst/images"
CACHE   <- "inst/extdata/asta_sez/cache"
options(fafbseg.flywire_dataset = "flywire_fafb_public")

SEZ <- list(
  CB0602  = "720575940632295751",  # SS32423 ∩ R65D05
  CB0239  = "720575940634548194"   # SS32423 ∩ VT019900
)
OTHERS <- list(
  SMP202_R = "720575940626412563",
  SMP202_L = "720575940628248326",
  SMP540_R = "720575940623972792",
  SMP540_L = "720575940616330011",
  CB2082_R = "720575940623952679",
  CB0108   = "720575940620093019"
)

# 1. Fetch / load meshes — cached.
sez_mesh_f    <- file.path(CACHE, "meshes_sez.rds")
others_mesh_f <- file.path(CACHE, "meshes_others.rds")

if (!file.exists(sez_mesh_f)) {
  cat("Fetching SEZ meshes (~20 s for 2 neurons) ...\n")
  m_sez <- fafbseg::read_cloudvolume_meshes(unname(unlist(SEZ)))
  saveRDS(m_sez, sez_mesh_f)
} else {
  m_sez <- readRDS(sez_mesh_f); cat("SEZ meshes cached\n")
}
names(m_sez) <- setNames(names(SEZ), unname(unlist(SEZ)))[names(m_sez)]

if (!file.exists(others_mesh_f)) {
  cat("Fetching non-SEZ meshes (~50 s for 6 neurons) ...\n")
  m_oth <- fafbseg::read_cloudvolume_meshes(unname(unlist(OTHERS)))
  saveRDS(m_oth, others_mesh_f)
} else {
  m_oth <- readRDS(others_mesh_f); cat("non-SEZ meshes cached\n")
}
names(m_oth) <- setNames(names(OTHERS), unname(unlist(OTHERS)))[names(m_oth)]

# 2. Brain + SEZ neuropil context.
brain    <- elmr::FAFB14.surf
sez_surf <- subset(elmr::FAFB14NP.surf,
                   intersect(c("GNG","SAD","FLA_R","FLA_L","AMMC_R","AMMC_L"),
                             elmr::FAFB14NP.surf$RegionList))

# geom_neuron.mesh3d uses geom_polygon (838k triangles) with high
# alpha — visually the polygons cancel each other out and the neuron
# disappears. For a denser, more MIP-like volumetric render, project
# the mesh vertices straight to 2D and lay them down as a low-alpha
# point cloud.
mesh_xy <- function(mesh) {
  v <- nat::xyzmatrix(mesh)
  data.frame(X = v[,1], Y = v[,2], Z = v[,3])
}
mk_panel <- function(mesh, title, neuron_col) {
  pts <- mesh_xy(mesh)
  ggplot() +
    geom_neuron(brain,    cols = "grey75", alpha = 0.30) +
    geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.40) +
    geom_point(data = pts, aes(x = X, y = Y),
               colour = neuron_col, alpha = 0.04, size = 0.15) +
    coord_fixed() + scale_y_reverse() +
    guides(colour = "none", fill = "none") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# --- SEZ pair (Step 8 figure)
cat("\nrender SEZ pair (meshes) ...\n")
p_cb0602 <- mk_panel(m_sez[["CB0602"]],
                     "CB0602 — SS32423 ∩ R65D05 (canonical AstA-GAL4)",
                     "magenta")
p_cb0239 <- mk_panel(m_sez[["CB0239"]],
                     "CB0239 — SS32423 ∩ VT019900 (LB11 hemilineage)",
                     "darkorange")
panel_sez <- (p_cb0602 | p_cb0239) +
  plot_annotation(
    title    = "Top-2 candidate FAFB-v783 neurons for the SEZ AstA cell",
    subtitle = "Frontal view, FAFB14 brain (grey) with SEZ neuropils (blue) — mesh renders"
  )
ggsave(file.path(OUT_DIR, "asta_sez_candidates_natggplot.png"),
       panel_sez, width = 11, height = 5.5, dpi = 300, bg = "white")
cat("wrote asta_sez_candidates_natggplot.png\n")

p_overlay <- ggplot() +
  geom_neuron(brain,    cols = "grey75", alpha = 0.30) +
  geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.40) +
  geom_point(data = mesh_xy(m_sez[["CB0602"]]), aes(x=X, y=Y),
             colour = "magenta",   alpha = 0.04, size = 0.15) +
  geom_point(data = mesh_xy(m_sez[["CB0239"]]), aes(x=X, y=Y),
             colour = "darkorange", alpha = 0.04, size = 0.15) +
  coord_fixed() + scale_y_reverse() +
  guides(colour = "none", fill = "none") +
  labs(title = "Top-2 candidates overlaid (CB0602 magenta, CB0239 orange)",
       subtitle = "FAFB14 brain (grey), SEZ neuropils (blue) — mesh renders") +
  theme_void(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(file.path(OUT_DIR, "asta_sez_candidates_overlay.png"),
       p_overlay, width = 7, height = 6, dpi = 300, bg = "white")
cat("wrote asta_sez_candidates_overlay.png\n")

# --- Non-SEZ panel (Step 10 figure)
cat("\nrender non-SEZ candidates (meshes) ...\n")
p1 <- mk_panel(m_oth[["SMP202_L"]], "SMP202 (L) — ACh, DCV pct 0.97",  "#08306b")
p2 <- mk_panel(m_oth[["SMP202_R"]], "SMP202 (R) — ACh, DCV pct 0.95",  "#08306b")
p3 <- mk_panel(m_oth[["SMP540_L"]], "SMP540 (L, SMPpd2) — glu, DCV pct 0.92", "#a6611a")
p4 <- mk_panel(m_oth[["SMP540_R"]], "SMP540 (R, SMPpd2) — glu, DCV pct 0.91", "#a6611a")
p5 <- mk_panel(m_oth[["CB2082_R"]], "CB2082 (SMPp&v1) — glu, DCV pct 0.93", "#018571")
p6 <- mk_panel(m_oth[["CB0108"]],   "CB0108 (LB19) — ACh, DCV pct 0.99 (IPS soma)", "#7b3294")

panel_oth <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title    = "Other SS32423-labelled candidates outside the SEZ",
    subtitle = "Top non-SEZ-soma SS32423 hits — mesh renders; bilateral SMP202 + SMP540 + CB2082, plus CB0108 at the SEZ/IPS boundary"
  )
ggsave(file.path(OUT_DIR, "asta_other_candidates_natggplot.png"),
       panel_oth, width = 13, height = 11, dpi = 300, bg = "white")
cat("wrote asta_other_candidates_natggplot.png\n")

p_oth_overlay <- ggplot() +
  geom_neuron(brain,    cols = "grey75", alpha = 0.30) +
  geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.40) +
  geom_point(data = mesh_xy(m_oth[["SMP202_L"]]), aes(x=X,y=Y), colour="#08306b", alpha=0.04, size=0.15) +
  geom_point(data = mesh_xy(m_oth[["SMP202_R"]]), aes(x=X,y=Y), colour="#08306b", alpha=0.04, size=0.15) +
  geom_point(data = mesh_xy(m_oth[["SMP540_L"]]), aes(x=X,y=Y), colour="#a6611a", alpha=0.04, size=0.15) +
  geom_point(data = mesh_xy(m_oth[["SMP540_R"]]), aes(x=X,y=Y), colour="#a6611a", alpha=0.04, size=0.15) +
  geom_point(data = mesh_xy(m_oth[["CB2082_R"]]), aes(x=X,y=Y), colour="#018571", alpha=0.04, size=0.15) +
  geom_point(data = mesh_xy(m_oth[["CB0108"]]),   aes(x=X,y=Y), colour="#7b3294", alpha=0.04, size=0.15) +
  coord_fixed() + scale_y_reverse() +
  guides(colour = "none", fill = "none") +
  labs(title = "Non-SEZ SS32423 candidates — all 6 overlaid",
       subtitle = "SMP202 (navy), SMP540 (brown), CB2082 (teal), CB0108 (purple)") +
  theme_void(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(file.path(OUT_DIR, "asta_other_candidates_overlay.png"),
       p_oth_overlay, width = 7, height = 6, dpi = 300, bg = "white")
cat("wrote asta_other_candidates_overlay.png\n")

# --- Re-render the depth-encoded MIP panel with meshes
cat("\nrender depth-encoded MIP panel (meshes) ...\n")
# For the depth-encoded MIP panel: skip nat.ggplot's 2-stop gradient
# and render the mesh polygons directly with the full
# PsychedelicRainBow ramp, keyed on Z (depth).
DEPTH_LUT <- c("#0033FF","#00CCFF","#00FF66","#CCFF00","#FFCC00",
               "#FF3300","#FF00CC")
mk_depth <- function(mesh, title) {
  pts <- mesh_xy(mesh)
  ggplot() +
    geom_neuron(brain, cols = "grey25", alpha = 0.18) +
    geom_point(data = pts, aes(x = X, y = Y, colour = Z),
               alpha = 0.06, size = 0.18) +
    scale_colour_gradientn(colours = DEPTH_LUT, guide = "none") +
    coord_fixed() + scale_y_reverse() +
    guides(colour = "none", fill = "none") +
    labs(title = title) +
    theme_void(base_size = 14) +
    theme(plot.background  = element_rect(fill = "black", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          plot.margin      = margin(6,6,6,6),
          plot.title       = element_text(hjust = 0.5, colour = "white",
                                          margin = margin(b = 4)),
          legend.position  = "none")
}
p_d602 <- mk_depth(m_sez[["CB0602"]], "")
p_d239 <- mk_depth(m_sez[["CB0239"]], "")

ggsave(file.path(OUT_DIR, "_tmp_cb0602.png"), p_d602,
       width = 14, height = 5.8, dpi = 200, bg = "black")
ggsave(file.path(OUT_DIR, "_tmp_cb0239.png"), p_d239,
       width = 14, height = 5.8, dpi = 200, bg = "black")

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
panel_mip <- magick::image_append(c(
  mk_hdr("SS32423  —  NeuronBridge brain MIP (v3_9_0, FlyLight Split-GAL4)"),
  ss_mip,
  mk_hdr("CB0602  —  SS32423 ∩ R65D05 (FAFB-v783 mesh, depth-encoded)", size = 30),
  cb0602_p,
  mk_hdr("CB0239  —  SS32423 ∩ VT019900 (FAFB-v783 mesh, depth-encoded)", size = 30),
  cb0239_p
), stack = TRUE)
magick::image_write(panel_mip,
                    file.path(OUT_DIR, "asta_sez_mip_panel.png"),
                    format = "png")
file.remove(file.path(OUT_DIR, c("_tmp_cb0602.png", "_tmp_cb0239.png")))
cat("wrote asta_sez_mip_panel.png\n")
