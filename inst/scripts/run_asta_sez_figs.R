#!/usr/bin/env Rscript
# Generate the AstA-SEZ figures from saved RDS, plus a 3D nat plot of the
# top candidates. Independent of the main pipeline so it's quick to iterate.

suppressMessages({
  library(arrow); library(dplyr); library(tidyr); library(ggplot2)
})
DATA_OUT <- file.path("inst", "extdata", "asta_sez")
OUT_DIR <- file.path("inst", "images")

ranked   <- readRDS(file.path(DATA_OUT, "asta_sez_ranked.rds"))
hits_top <- readRDS(file.path(DATA_OUT, "asta_sez_hits_top.rds"))
panel <- c("SS32423","R65D05","R65D06","R65D07","R65E01","VT019900")

# Re-derive a clean rank including the rank_score I dropped earlier.
ranked <- ranked %>% mutate(
  rank_score = (n_lines * 1.0)
             + (best_score / 50000)
             + ifelse(in_ss32423, 0.5, 0)
             + ifelse(coalesce(dcv_rich, FALSE), 0.5, 0)
             + ifelse(in_R65D05, 0.5, 0)
)
ranked <- ranked %>% arrange(desc(in_ss32423), desc(n_lines), desc(rank_score))

cat("\n========== TOP 25 SEZ-OK CANDIDATES (full meta) ==========\n")
print(ranked %>% head(25) %>%
        select(root_783, in_ss32423, in_R65D05, n_lines, lines, best_score,
               soma_zone, frac_sez, cell_type, hemilineage, super_class,
               nt, dcv_pct, dcv_rich, np_verified, rank_score),
      n = 25, width = Inf)

ss <- ranked %>% filter(in_ss32423) %>%
        arrange(desc(in_R65D05), desc(rank_score))
cat("\n\n========== ALL SS32423 HITS (sort by R65D05 + rank) ==========\n")
print(ss %>% select(root_783, in_R65D05, n_lines, lines, best_score,
                    soma_zone, cell_type, hemilineage, super_class, nt,
                    dcv_pct, dcv_rich, np_verified, rank_score),
      n = nrow(ss), width = Inf)

# --- ggplot summary
top15 <- ranked %>% head(15) %>%
  mutate(label = ifelse(!is.na(cell_type) & cell_type != "",
                        paste0(cell_type, "\n", substr(root_783, 11, 18)),
                        substr(root_783, 11, 18)),
         soma_zone = factor(soma_zone,
                            levels = c("SEZ","SEZ_borderline","unknown","non_SEZ")))

p_match <- ggplot(top15,
                  aes(x = n_lines, y = best_score,
                      colour = soma_zone,
                      size = pmax(coalesce(dcv_pct, 0), 0.05),
                      label = label)) +
  geom_point() +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf) +
  scale_colour_manual(values = c(SEZ="#08306b", SEZ_borderline="#6baed6",
                                 unknown="grey60", non_SEZ="#a50f15"),
                      drop = FALSE) +
  scale_size_continuous(range = c(2,7), limits = c(0,1),
                        labels = scales::percent, name = "DCV percentile") +
  labs(x = "n driver lines hitting this neuron",
       y = "best NeuronBridge score",
       colour = "soma zone",
       title = "Candidate FAFB-783 neurons for the SEZ AstA cell",
       subtitle = "AstA-locus driver consensus × NeuronBridge × DCV soma density (top 15)") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT_DIR, "asta_sez_match.png"), p_match,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("wrote asta_sez_match.png\n")

p_consensus <- top15 %>%
  ggplot(aes(x = reorder(label, n_lines), y = n_lines,
             fill = coalesce(dcv_rich, FALSE))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE`="#08306b", `FALSE`="grey70"),
                    labels = c(`TRUE`="DCV-rich (>=p90)",
                               `FALSE`="below p90"),
                    name = NULL) +
  labs(x = NULL, y = "n driver lines",
       title = "Cross-line consensus (top 15)") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT_DIR, "asta_sez_consensus.png"), p_consensus,
       width = 8, height = 6, dpi = 300, bg = "white")
cat("wrote asta_sez_consensus.png\n")

# Heatmap of which of our 6 panel lines hit each top-25 candidate
top25_ids <- head(ranked$root_783, 25)
mat_long <- hits_top %>%
  filter(root_783 %in% top25_ids) %>%
  mutate(in_hit = 1) %>%
  select(query_line, root_783, in_hit) %>%
  complete(query_line = panel,
           root_783 = top25_ids,
           fill = list(in_hit = 0)) %>%
  left_join(ranked %>% select(root_783, cell_type), by = "root_783") %>%
  mutate(query_line = factor(query_line, levels = panel),
         label = ifelse(!is.na(cell_type) & cell_type != "",
                        paste0(cell_type, "·", substr(root_783, 14, 18)),
                        substr(root_783, 11, 18)),
         label = factor(label,
                        levels = unique(setNames(label, root_783)[top25_ids])))
p_heat <- ggplot(mat_long,
                 aes(x = query_line, y = label, fill = factor(in_hit))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(`0`="grey90", `1`="#08306b"), guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = "candidate FlyWire neuron",
       title = "Driver-panel hit matrix — top 25 candidates",
       subtitle = "Filled = line hits this neuron in NeuronBridge v3_9_0 (cap 10 MIPs/line)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8))
ggsave(file.path(OUT_DIR, "asta_sez_heatmap.png"), p_heat,
       width = 7.5, height = 9, dpi = 300, bg = "white")
cat("wrote asta_sez_heatmap.png\n")
