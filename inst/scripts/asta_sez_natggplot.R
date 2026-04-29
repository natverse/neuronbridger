#!/usr/bin/env Rscript
# Render the top-2 candidate FAFB-v783 neurons as a 2D nat.ggplot
# figure on the FAFB14 brain with SEZ neuropils outlined.
# Mirrors the dcv-repo pattern (R/visualise/fafb_flange.R,
# fig_2_dcv_predictions_fafb.Rmd).

suppressMessages({
  library(nat); library(nat.ggplot); library(nat.flybrains); library(elmr)
  library(nat.templatebrains); library(fafbseg)
  library(ggplot2); library(dplyr); library(patchwork)
})

OUT_DIR <- "inst/images"
CACHE   <- "inst/extdata/asta_sez/cache"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CB0602  <- "720575940632295751"  # SS32423 ∩ R65D05
CB0239  <- "720575940634548194"  # SS32423 ∩ VT019900

# 1. L2 skeletons (FlyWire raw → nm in FAFB14 space).
sk_path <- file.path(CACHE, "l2skel_top2.rds")
if (!file.exists(sk_path)) {
  sk <- read_l2skel(c(CB0602, CB0239),
                    datastack_name = "flywire_fafb_public")
  saveRDS(sk, sk_path)
} else {
  sk <- readRDS(sk_path)
}
names(sk) <- c("CB0602", "CB0239")[match(names(sk), c(CB0602, CB0239))]

# 2. Brain meshes — FAFB14 native (FlyWire uses FAFB14 nm), and SEZ
# neuropils from JFRC2NP transformed in.
brain  <- elmr::FAFB14.surf
np_all <- elmr::FAFB14NP.surf
SEZ <- intersect(c("GNG","SAD","FLA_R","FLA_L","AMMC_R","AMMC_L"),
                 np_all$RegionList)
sez_surf <- subset(np_all, SEZ)

# 3. ggplot panels — frontal view (XY, flip Y so dorsal is up).
mk_panel <- function(neuron, title, neuron_col) {
  ggplot() +
    geom_neuron(brain, cols = "grey70", alpha = 0.35) +
    geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.45) +
    geom_neuron(neuron, cols = neuron_col, alpha = 0.95, lwd = 0.4) +
    coord_fixed() +
    scale_y_reverse() +
    guides(colour = "none", fill = "none") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

p1 <- mk_panel(sk[["CB0602"]],
               "CB0602 — SS32423 ∩ R65D05 (canonical AstA-GAL4)",
               "magenta")
p2 <- mk_panel(sk[["CB0239"]],
               "CB0239 — SS32423 ∩ VT019900 (LB11 hemilineage)",
               "darkorange")

panel <- (p1 | p2) +
  plot_annotation(
    title = "Top-2 candidate FAFB-v783 neurons for the SEZ AstA cell",
    subtitle = "Frontal view, FAFB14 brain (grey) with SEZ neuropils (blue)"
  )
ggsave(file.path(OUT_DIR, "asta_sez_candidates_natggplot.png"),
       panel, width = 11, height = 5.5, dpi = 300, bg = "white")
cat("wrote asta_sez_candidates_natggplot.png\n")

# 4. Combined panel: both neurons overlaid in same view (so they can be
# directly compared with the SS32423 NB MIP).
p_both <- ggplot() +
  geom_neuron(brain, cols = "grey70", alpha = 0.35) +
  geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.45) +
  geom_neuron(sk[["CB0602"]], cols = "magenta", alpha = 0.95, lwd = 0.4) +
  geom_neuron(sk[["CB0239"]], cols = "darkorange", alpha = 0.95, lwd = 0.4) +
  coord_fixed() +
  scale_y_reverse() +
  guides(colour = "none", fill = "none") +
  labs(title = "Top-2 candidates overlaid (CB0602 magenta, CB0239 orange)",
       subtitle = "FAFB14 brain (grey), SEZ neuropils (blue)") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(file.path(OUT_DIR, "asta_sez_candidates_overlay.png"),
       p_both, width = 7, height = 6, dpi = 300, bg = "white")
cat("wrote asta_sez_candidates_overlay.png\n")
