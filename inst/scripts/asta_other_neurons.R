#!/usr/bin/env Rscript
# Beyond the SEZ: render the other AstA-cell candidates SS32423 labels.
# Inspecting the SS32423 MIP montage (`asta_sez_SS32423_montage.png`),
# row 1 cols 2–4 show bright cells outside the SEZ — and the 16
# non-SEZ-soma SS32423 hits in `ranked_full.csv` cluster strongly in
# the SMP / SLP region with several bilateral cell-type recurrences.
# These are good candidates for the AstA cells outside the SEZ
# described by Hentze 2015 / Pool 2014 (PI / SMP / pars lateralis).
#
# We render the strongest non-SEZ candidates as a nat.ggplot panel and
# stack them under the SEZ candidates for direct comparison.

suppressMessages({
  library(devtools); load_all(".")
  library(nat); library(nat.ggplot); library(elmr); library(fafbseg)
  library(dplyr); library(ggplot2); library(patchwork); library(magick)
})

OUT_DIR <- "inst/images"
CACHE   <- "inst/extdata/asta_sez/cache"

# Pick the strongest non-SEZ SS32423 candidates from ranked_full.csv.
# Prefer DCV-rich + bilateral recurrence (same cell_type appearing for
# both sides). These were the bilateral pairs identified at run time:
#   SMP202   (×2, ACh,    DCV pct 0.95/0.97)
#   SMP540   (×3, gluta., DCV pct 0.91)
#   CB2082   (×2, gluta., DCV pct 0.93/0.86)
#   CB0108   (×1, ACh,    DCV pct 0.99)  — IPS-adjacent SEZ cell
OTHERS <- list(
  SMP202_R = "720575940626412563",  # right SMP202, ACh, DCV-rich, dcv_pct 0.95
  SMP202_L = "720575940628248326",  # left  SMP202, ACh, DCV-rich, dcv_pct 0.97
  SMP540_R = "720575940623972792",  # right SMP540, glu, DCV-rich, dcv_pct 0.91
  SMP540_L = "720575940616330011",  # left  SMP540, glu, DCV-rich, dcv_pct 0.92
  CB2082_R = "720575940623952679",  # right CB2082, glu, DCV-rich, dcv_pct 0.93
  CB0108   = "720575940620093019"   # CB0108,   ACh, DCV-rich (0.99), IPS soma
)

# 1. Fetch L2 skeletons (cached).
sk_path <- file.path(CACHE, "l2skel_others.rds")
if (!file.exists(sk_path)) {
  sk <- read_l2skel(unname(unlist(OTHERS)),
                    datastack_name = "flywire_fafb_public")
  saveRDS(sk, sk_path)
} else {
  sk <- readRDS(sk_path)
}
# Map FlyWire root → readable name for the panels.
nm_lookup <- setNames(names(OTHERS), unname(unlist(OTHERS)))
names(sk)  <- nm_lookup[names(sk)]

# 2. Brain mesh and SEZ neuropils for context.
brain  <- elmr::FAFB14.surf
sez_surf <- subset(elmr::FAFB14NP.surf,
                   intersect(c("GNG","SAD","FLA_R","FLA_L","AMMC_R","AMMC_L"),
                             elmr::FAFB14NP.surf$RegionList))

mk_panel <- function(neuron, title, neuron_col) {
  ggplot() +
    geom_neuron(brain,    cols = "grey70", alpha = 0.30) +
    geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.35) +
    geom_neuron(neuron,   cols = neuron_col, alpha = 0.95, lwd = 0.4) +
    coord_fixed() + scale_y_reverse() +
    guides(colour = "none", fill = "none") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "none")
}

p1 <- mk_panel(sk[["SMP202_L"]], "SMP202 (L) — ACh, DCV pct 0.97",  "#08306b")
p2 <- mk_panel(sk[["SMP202_R"]], "SMP202 (R) — ACh, DCV pct 0.95",  "#08306b")
p3 <- mk_panel(sk[["SMP540_L"]], "SMP540 (L, SMPpd2) — glu, DCV pct 0.92", "#a6611a")
p4 <- mk_panel(sk[["SMP540_R"]], "SMP540 (R, SMPpd2) — glu, DCV pct 0.91", "#a6611a")
p5 <- mk_panel(sk[["CB2082_R"]], "CB2082 (SMPp&v1) — glu, DCV pct 0.93", "#018571")
p6 <- mk_panel(sk[["CB0108"]],   "CB0108 (LB19) — ACh, DCV pct 0.99 (IPS soma)", "#7b3294")

panel <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title    = "Other SS32423-labelled candidates outside the SEZ",
    subtitle = "Top non-SEZ-soma SS32423 hits — bilateral SMP202 + SMP540 + CB2082, plus CB0108 at the SEZ/IPS boundary"
  )
ggsave(file.path(OUT_DIR, "asta_other_candidates_natggplot.png"), panel,
       width = 13, height = 11, dpi = 300, bg = "white")
cat("wrote asta_other_candidates_natggplot.png\n")

# 3. Combined overlay: all 5 candidates on one brain — gives a sense of
# the *spatial spread* of the non-SEZ AstA-candidate set.
DEPTH_LUT <- c("#0033FF","#00CCFF","#00FF66","#CCFF00","#FFCC00",
               "#FF3300","#FF00CC")
p_overlay <- ggplot() +
  geom_neuron(brain,    cols = "grey70", alpha = 0.30) +
  geom_neuron(sez_surf, cols = "#3182bd", alpha = 0.35) +
  geom_neuron(sk[["SMP202_L"]], cols = "#08306b", alpha = 0.85, lwd = 0.4) +
  geom_neuron(sk[["SMP202_R"]], cols = "#08306b", alpha = 0.85, lwd = 0.4) +
  geom_neuron(sk[["SMP540_L"]], cols = "#a6611a", alpha = 0.85, lwd = 0.4) +
  geom_neuron(sk[["SMP540_R"]], cols = "#a6611a", alpha = 0.85, lwd = 0.4) +
  geom_neuron(sk[["CB2082_R"]], cols = "#018571", alpha = 0.85, lwd = 0.4) +
  geom_neuron(sk[["CB0108"]],   cols = "#7b3294", alpha = 0.85, lwd = 0.4) +
  coord_fixed() + scale_y_reverse() +
  guides(colour = "none", fill = "none") +
  labs(title = "Non-SEZ SS32423 candidates — all 5 overlaid",
       subtitle = "SMP202 (navy), SMP540 (brown), CB2082 (teal), CB0108 (purple)") +
  theme_void(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(file.path(OUT_DIR, "asta_other_candidates_overlay.png"),
       p_overlay, width = 7, height = 6, dpi = 300, bg = "white")
cat("wrote asta_other_candidates_overlay.png\n")
