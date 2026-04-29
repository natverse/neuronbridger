#!/usr/bin/env Rscript
# Reproducer for vignettes/asta_sez.Rmd
# Caches NB hits per line in inst/extdata/asta_sez/cache/ so re-runs are
# fast. Writes ranked CSV + ggplot summary PNGs into inst/images/.

suppressMessages({
  library(devtools); load_all(".")
  library(arrow); library(dplyr); library(tidyr); library(ggplot2)
})
NB_VERSION <- "v3_9_0"
DCV_DIR <- "/Users/GD/LMBD/Papers/dcv/data/sjcabs/fafb"
OUT_DIR <- file.path("inst", "images")
DATA_OUT <- file.path("inst", "extdata", "asta_sez")
CACHE <- file.path(DATA_OUT, "cache")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

driver_panel <- tibble::tribble(
  ~line,        ~role,
  "SS32423",    "primary lead — aMulet split (Sterne et al. 2021)",
  "R65D05",     "Pfeiffer 'AstA-GAL4' (Hergarden 2012)",
  "R65D06",     "DBD half of SS32423",
  "R65D04",     "neighbouring AstA-locus tile",
  "R65D07",     "neighbouring AstA-locus tile",
  "R65E01",     "neighbouring AstA-locus tile",
  "VT019900",   "AD half of SS32423"
)

# Reduce hit-row payload to a stable narrow schema before binding,
# so per-MIP variable annotation columns don't blow up bind_rows.
KEEP_COLS <- c("publishedName","libraryName","normalizedScore",
               "matchingPixels","alignmentSpace","anatomicalArea","gender")

cat("Step 1: confirm panel in NB", NB_VERSION, "\n")
panel_info <- list()
for (ln in driver_panel$line) {
  out <- try(neuronbridge_info(ln, dataset = "by_line", version = NB_VERSION), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out) || !nrow(out)) {
    cat("  -", ln, ": NOT IN", NB_VERSION, "\n"); next
  }
  brain <- out[out$anatomicalArea == "Brain", c("publishedName","nb.id","libraryName","alignmentSpace","anatomicalArea","gender"), drop = FALSE]
  # Cap MIPs per line — broad GAL4 lines have many MIPs and each
  # neuronbridge_hits call uses plyr::rbind.fill which is quadratic in
  # the number of result rows. Multiple MIPs of one line are
  # near-redundant (same GAL4, slightly different staining) and we
  # already best-MIP-per-line downstream, so 5 is plenty.
  MIPS_PER_LINE_CAP <- 10L
  if (nrow(brain) > MIPS_PER_LINE_CAP) brain <- head(brain, MIPS_PER_LINE_CAP)
  cat(sprintf("  - %-10s %d MIPs total (%d brain, kept %d)\n",
              ln, nrow(out), sum(out$anatomicalArea == "Brain"), nrow(brain)))
  brain$line <- ln
  panel_info[[ln]] <- brain
}
panel_info <- bind_rows(panel_info)
cat("  total brain MIPs across panel:", nrow(panel_info), "\n")

cat("\nStep 2: pull FlyWire-FAFB-v783 hits per brain MIP (cached per line)\n")
hits_list <- list()
for (ln in unique(panel_info$line)) {
  cache_f <- file.path(CACHE, paste0("hits_", ln, ".rds"))
  if (file.exists(cache_f)) {
    cat("  -", ln, ": cached\n")
    hits_list[[ln]] <- readRDS(cache_f)
    next
  }
  ln_mips <- panel_info[panel_info$line == ln, ]
  rows <- list()
  for (i in seq_len(nrow(ln_mips))) {
    nbid <- ln_mips$nb.id[i]
    h <- try(neuronbridge_hits(nbid, version = NB_VERSION), silent = TRUE)
    if (inherits(h, "try-error") || is.null(h) || !nrow(h)) next
    h <- h[grepl("FlyWire", h$libraryName), , drop = FALSE]
    if (!nrow(h)) next
    # narrow to stable cols only
    keep <- intersect(KEEP_COLS, colnames(h))
    h <- as.data.frame(h[, keep, drop = FALSE])
    h$normalizedScore <- as.numeric(h$normalizedScore)
    h$matchingPixels  <- as.numeric(h$matchingPixels)
    h$query_line <- ln
    h$query_mip  <- nbid
    rows[[length(rows)+1L]] <- h
  }
  ln_hits <- if (length(rows)) bind_rows(rows) else tibble()
  saveRDS(ln_hits, cache_f)
  cat(sprintf("  - %-10s %4d FlyWire hits (across %d brain MIPs)\n",
              ln, nrow(ln_hits), nrow(ln_mips)))
  hits_list[[ln]] <- ln_hits
}
hits_raw <- bind_rows(hits_list)
hits_raw$root_783 <- sub("^flywire_fafb:v783:", "", hits_raw$publishedName)
cat("  total FlyWire hits across panel:", nrow(hits_raw), "\n")

# Best-MIP-per-line-per-neuron
hits_best <- hits_raw %>%
  group_by(query_line, root_783) %>%
  slice_max(normalizedScore, n = 1, with_ties = FALSE) %>%
  ungroup()

# Per-line elbow cap (don't cap too low — we want the union to be wide)
top_with_elbow <- function(df, cap = 25, drop_frac = 0.75, floor_frac = 0.20) {
  df <- df[order(df$normalizedScore, decreasing = TRUE), , drop = FALSE]
  if (!nrow(df)) return(df)
  keep <- TRUE
  if (nrow(df) > 1) {
    ratios <- df$normalizedScore[-1] / df$normalizedScore[-nrow(df)]
    keep <- c(TRUE,
              cumprod(ratios >= drop_frac) == 1 &
              df$normalizedScore[-1] / df$normalizedScore[1] >= floor_frac)
  }
  head(df[keep, , drop = FALSE], cap)
}
hits_top <- hits_best %>% group_by(query_line) %>%
  group_modify(~ top_with_elbow(.x, cap = 25)) %>% ungroup()
cat("  hits_top rows after elbow:", nrow(hits_top),
    " unique neurons:", n_distinct(hits_top$root_783), "\n")

cat("\nStep 3: soma-neuropil tokenisation from soma-DCV detection feather\n")
SEZ_INNER <- c("GNG","SAD","AMMC_L","AMMC_R","FLA_L","FLA_R")
SEZ_TOK <- c(SEZ_INNER, paste0("outside_", SEZ_INNER))
soma_class_f <- file.path(CACHE, "soma_class.rds")
if (file.exists(soma_class_f)) {
  cat("  cached\n"); soma_class <- readRDS(soma_class_f)
} else {
  soma_dcv <- read_feather(file.path(DCV_DIR, "fafb_783_soma_dcv_detection.feather"),
                           col_select = c("root_783","neuropil"))
  sez_tok <- function(np, set) vapply(strsplit(np, ",", fixed=TRUE),
                                      function(x) any(x %in% set), logical(1))
  soma_dcv <- soma_dcv %>%
    mutate(any_sez   = sez_tok(neuropil, SEZ_TOK),
           inner_sez = sez_tok(neuropil, SEZ_INNER))
  soma_class <- soma_dcv %>% group_by(root_783) %>%
    summarise(n_dcv = n(),
              frac_sez = mean(any_sez),
              frac_sez_inner = mean(inner_sez),
              top_token = names(sort(table(neuropil), decreasing=TRUE))[1],
              .groups="drop") %>%
    mutate(soma_zone = case_when(frac_sez >= 0.5 ~ "SEZ",
                                 frac_sez >= 0.1 ~ "SEZ_borderline",
                                 TRUE ~ "non_SEZ"))
  saveRDS(soma_class, soma_class_f)
}
cat("  soma_zone counts:\n"); print(table(soma_class$soma_zone))

cat("\nStep 4: DCV-density percentiles (central-brain baseline)\n")
meta <- read_feather(file.path(DCV_DIR, "fafb_783_meta.feather")) %>%
  mutate(root_783 = as.character(fafb_783_id))
cb <- meta %>% filter(region == "central_brain", !is.na(soma_dcv_density))
dcv_thr <- quantile(cb$soma_dcv_density, probs = 0.90, na.rm = TRUE)
ecdf_cb <- ecdf(cb$soma_dcv_density)
meta <- meta %>% mutate(dcv_pct  = ecdf_cb(soma_dcv_density),
                        dcv_rich = soma_dcv_density >= dcv_thr)
cat("  p90 threshold on central-brain soma_dcv_density:",
    sprintf("%.3f", dcv_thr), "\n")

cat("\nStep 5: rank\n")
cand <- hits_top %>%
  left_join(soma_class, by = "root_783") %>%
  left_join(meta %>% select(root_783, cell_class, cell_sub_class, cell_type,
                            super_class, hemilineage, side, flow,
                            neurotransmitter_predicted, neuropeptide_verified,
                            soma_dcv_density, dcv_pct, dcv_rich),
            by = "root_783")

# Use base-R aggregation to dodge dplyr::first() quirks on heterogeneous cols.
ranked <- cand %>%
  mutate(in_ss32423 = query_line == "SS32423",
         in_R65D05  = query_line == "R65D05") %>%
  group_by(root_783) %>%
  summarise(
    n_lines       = n_distinct(query_line),
    in_ss32423    = any(in_ss32423),
    in_R65D05     = any(in_R65D05),
    lines         = paste(sort(unique(query_line)), collapse = ","),
    best_score    = max(normalizedScore, na.rm = TRUE),
    score_sum     = sum(pmin(normalizedScore, 50000), na.rm = TRUE),
    soma_zone     = soma_zone[1],
    frac_sez      = frac_sez[1],
    top_token     = top_token[1],
    soma_dcv_density = soma_dcv_density[1],
    dcv_pct       = dcv_pct[1],
    dcv_rich      = dcv_rich[1],
    cell_class    = cell_class[1],
    cell_type     = cell_type[1],
    hemilineage   = hemilineage[1],
    super_class   = super_class[1],
    nt            = neurotransmitter_predicted[1],
    np_verified   = neuropeptide_verified[1],
    .groups = "drop"
  ) %>%
  mutate(soma_zone = ifelse(is.na(soma_zone), "unknown", soma_zone),
         sez_ok    = soma_zone %in% c("SEZ","SEZ_borderline","unknown"),
         dcv_rich  = ifelse(is.na(dcv_rich), FALSE, dcv_rich),
         rank_score = (n_lines * 1.0)
                    + (best_score / 50000)
                    + ifelse(in_ss32423, 0.5, 0)
                    + ifelse(dcv_rich, 0.5, 0)
                    + ifelse(in_R65D05, 0.5, 0)) %>%
  arrange(desc(in_ss32423), desc(n_lines), desc(rank_score))

cat("  total candidate FlyWire neurons:", nrow(ranked), "\n")
cat("  in SS32423:", sum(ranked$in_ss32423), "\n")
cat("  also in R65D05:", sum(ranked$in_ss32423 & ranked$in_R65D05), "\n")

ss <- ranked %>% filter(in_ss32423)
cat(sprintf("\nSANITY: SS32423 candidates retained: %d  (SEZ=%d  borderline=%d  unknown=%d  non_SEZ=%d)\n",
    nrow(ss),
    sum(ss$soma_zone == "SEZ"),
    sum(ss$soma_zone == "SEZ_borderline"),
    sum(ss$soma_zone == "unknown"),
    sum(ss$soma_zone == "non_SEZ")))

ranked_sez <- ranked %>% filter(sez_ok)
write.csv(ranked,     file.path(DATA_OUT, "asta_sez_ranked_full.csv"),  row.names = FALSE)
write.csv(head(ranked_sez, 30), file.path(DATA_OUT, "asta_sez_ranked_top30.csv"), row.names = FALSE)
saveRDS(ranked_sez, file.path(DATA_OUT, "asta_sez_ranked.rds"))
saveRDS(hits_top,   file.path(DATA_OUT, "asta_sez_hits_top.rds"))

cat("\n========== TOP 25 SEZ-OK CANDIDATES ==========\n")
print(head(ranked_sez %>% select(root_783, in_ss32423, in_R65D05, n_lines, lines,
                                 best_score, soma_zone, frac_sez,
                                 cell_type, cell_class, hemilineage, super_class, nt,
                                 dcv_pct, dcv_rich, np_verified), 25),
      n = 25, width = Inf)

cat("\n\n========== ALL SS32423 HITS (sorted by R65D05+rank) ==========\n")
ss_sort <- ss %>% arrange(desc(in_R65D05), desc(rank_score))
print(ss_sort %>% select(root_783, in_R65D05, n_lines, lines, best_score,
                         soma_zone, frac_sez, top_token,
                         cell_type, cell_class, hemilineage, super_class, nt,
                         dcv_pct, dcv_rich, np_verified, rank_score),
      n = nrow(ss_sort), width = Inf)

cat("\n\nStep 6: ggplot summary figures\n")
top15 <- head(ranked_sez, 15) %>%
  mutate(label = ifelse(!is.na(cell_type) & cell_type != "",
                        paste0(cell_type, "\n", substr(root_783, 11, 18)),
                        substr(root_783, 11, 18)),
         soma_zone = factor(soma_zone, levels = c("SEZ","SEZ_borderline","unknown","non_SEZ")))

p_match <- ggplot(top15,
                  aes(x = n_lines, y = best_score,
                      colour = soma_zone, size = pmax(dcv_pct, 0.05, na.rm = TRUE),
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

p_consensus <- top15 %>%
  ggplot(aes(x = reorder(label, n_lines), y = n_lines,
             fill = dcv_rich)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE`="#08306b", `FALSE`="grey70"),
                    labels = c(`TRUE`="DCV-rich (>=p90)", `FALSE`="below p90"),
                    name = NULL) +
  labs(x = NULL, y = "n driver lines",
       title = "Cross-line consensus (top 15)") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT_DIR, "asta_sez_consensus.png"), p_consensus,
       width = 8, height = 6, dpi = 300, bg = "white")

# Heatmap of which lines hit the top 20 candidates
mat_long <- hits_top %>% mutate(in_hit = 1) %>%
  filter(root_783 %in% head(ranked_sez$root_783, 20)) %>%
  select(query_line, root_783, in_hit) %>%
  complete(query_line = unique(driver_panel$line),
           root_783 = head(ranked_sez$root_783, 20),
           fill = list(in_hit = 0)) %>%
  mutate(query_line = factor(query_line, levels = driver_panel$line),
         root_short = substr(root_783, 11, 18),
         root_short = factor(root_short,
                             levels = unique(substr(head(ranked_sez$root_783, 20), 11, 18))))
p_heat <- ggplot(mat_long, aes(x = query_line, y = root_short, fill = factor(in_hit))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(`0`="grey90", `1`="#08306b"), guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = "candidate FlyWire neuron (last 8 digits)",
       title = "Driver-panel hit matrix — top 20 candidates",
       subtitle = "Filled = line hits this neuron in NeuronBridge v3_9_0") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "asta_sez_heatmap.png"), p_heat,
       width = 7, height = 7, dpi = 300, bg = "white")

cat("\nAll outputs written to", OUT_DIR, "and", DATA_OUT, "\n")
