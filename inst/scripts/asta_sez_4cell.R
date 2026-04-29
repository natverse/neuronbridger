#!/usr/bin/env Rscript
# Render the canonical 4-cell SEZ AstA pattern from Hergarden 2012
# (PMC3309792) using FlyWire-FAFB-v783 meshes:
#   CB0602  (left + right) — ventral SEZ pair, putative_primary, ACh
#   CB0108  (left + right) — dorsal/IPS  pair, LB19,             ACh
# All 4 are SS32423 hits in the raw caches; CB0602_R + CB0108_L fell
# below the original elbow-cap threshold but appear in the per-line
# RDS caches.

suppressMessages({
  library(devtools); load_all(".")
  library(nat); library(nat.ggplot); library(elmr); library(fafbseg)
  library(ggplot2); library(patchwork); library(magick)
})
options(fafbseg.flywire_dataset = "flywire_fafb_public")

OUT_DIR <- "inst/images"
CACHE   <- "inst/extdata/asta_sez/cache"

FOUR <- list(
  CB0602_L = "720575940632295751",  # SS32423 ∩ R65D05 (top hit)
  CB0602_R = "720575940640469848",  # SS32423 + R65D05 + R65D06 caches (sub-elbow)
  CB0108_L = "720575940624028855",  # SS32423 + R65D06 caches
  CB0108_R = "720575940620093019"   # SS32423 (DCV pct 0.99); IPS soma
)

# Fetch / cache meshes.
mesh_f <- file.path(CACHE, "meshes_4cell.rds")
if (!file.exists(mesh_f)) {
  cat("Fetching 4 FlyWire meshes (~40 s) ...\n")
  m <- fafbseg::read_cloudvolume_meshes(unname(unlist(FOUR)))
  saveRDS(m, mesh_f)
} else {
  m <- readRDS(mesh_f); cat("4-cell meshes cached\n")
}
names(m) <- setNames(names(FOUR), unname(unlist(FOUR)))[names(m)]

brain    <- elmr::FAFB14.surf
sez_surf <- subset(elmr::FAFB14NP.surf,
                   intersect(c("GNG","SAD","FLA_R","FLA_L","AMMC_R","AMMC_L",
                               "PRW","IPS_L","IPS_R","SPS_L","SPS_R",
                               "VES_L","WED_L","WED_R"),
                             elmr::FAFB14NP.surf$RegionList))

mesh_xy <- function(mesh) {
  v <- nat::xyzmatrix(mesh)
  data.frame(X = v[,1], Y = v[,2], Z = v[,3])
}
mk_panel <- function(mesh, title, neuron_col) {
  ggplot() +
    geom_neuron(brain,    cols = "grey75", alpha = 0.30) +
    geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.40) +
    geom_point(data = mesh_xy(mesh), aes(x = X, y = Y),
               colour = neuron_col, alpha = 0.04, size = 0.15) +
    coord_fixed() + scale_y_reverse() +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# 2x2 panel — bilateral pairs side by side, dorsal pair on top
p_l108 <- mk_panel(m[["CB0108_L"]], "CB0108 (L) — dorsal SEZ/IPS pair",     "#7b3294")
p_r108 <- mk_panel(m[["CB0108_R"]], "CB0108 (R) — dorsal SEZ/IPS pair",     "#7b3294")
p_l602 <- mk_panel(m[["CB0602_L"]], "CB0602 (L) — ventral SEZ pair",        "magenta")
p_r602 <- mk_panel(m[["CB0602_R"]], "CB0602 (R) — ventral SEZ pair",        "magenta")
panel <- (p_l108 | p_r108) / (p_l602 | p_r602) +
  plot_annotation(
    title    = "The 4 canonical SEZ AstA cells (Hergarden 2012) — FAFB-v783",
    subtitle = "Two bilateral pairs: CB0108 (dorsal/IPS, LB19) above CB0602 (ventral SEZ, putative_primary)"
  )
ggsave(file.path(OUT_DIR, "asta_sez_4cell_panel.png"),
       panel, width = 11, height = 11, dpi = 300, bg = "white")
cat("wrote asta_sez_4cell_panel.png\n")

# Combined overlay — all 4 cells on the same brain, colour-coded by
# pair so readers see the bilateral symmetry + dorsal/ventral split.
p_overlay <- ggplot() +
  geom_neuron(brain,    cols = "grey75", alpha = 0.25) +
  geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.35) +
  geom_point(data = mesh_xy(m[["CB0602_L"]]), aes(x=X,y=Y),
             colour = "magenta",  alpha = 0.04, size = 0.18) +
  geom_point(data = mesh_xy(m[["CB0602_R"]]), aes(x=X,y=Y),
             colour = "magenta",  alpha = 0.04, size = 0.18) +
  geom_point(data = mesh_xy(m[["CB0108_L"]]), aes(x=X,y=Y),
             colour = "#7b3294",  alpha = 0.04, size = 0.18) +
  geom_point(data = mesh_xy(m[["CB0108_R"]]), aes(x=X,y=Y),
             colour = "#7b3294",  alpha = 0.04, size = 0.18) +
  coord_fixed() + scale_y_reverse() +
  labs(title = "All 4 SEZ AstA cells overlaid on the FAFB14 brain",
       subtitle = "Magenta = CB0602 (ventral pair); Purple = CB0108 (dorsal/IPS pair); SEZ neuropils blue") +
  theme_void(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(file.path(OUT_DIR, "asta_sez_4cell_overlay.png"),
       p_overlay, width = 7, height = 6, dpi = 300, bg = "white")
cat("wrote asta_sez_4cell_overlay.png\n")
