#!/usr/bin/env Rscript
# Reproducer for vignettes/abdominal_peripheral_targets.Rmd:
# pulls MANC abdominal motor + endocrine neurons, fetches NB v3_9_0 hits,
# expands SS/IS splits to AD+DBD hemidrivers, runs KDRC lookups, and
# renders the figures embedded in the vignette into inst/images/.
#
# KDRC chromote queries are slow (~1 hour for 500 lines). Pass
# `--skip-kdrc` to skip them and use a small example KDRC subset for
# demonstration plotting.

suppressMessages({
  library(devtools); load_all(".")
  library(malevnc); library(neuprintr); library(nat)
  library(dplyr); library(tidyr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
SKIP_KDRC <- "--skip-kdrc" %in% args
NB_VERSION <- "v3_9_0"
OUT_DIR <- "inst/images"
DATA_OUT <- "inst/extdata/abdominal"
CACHE   <- file.path(DATA_OUT, "cache")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE,   showWarnings = FALSE, recursive = TRUE)

cat("Step 1: pull MANC abdominal targets\n")
targets_f <- file.path(CACHE, "manc_targets.rds")
if (file.exists(targets_f)) {
  targets <- readRDS(targets_f)
  cat("  cached: nrow=", nrow(targets), "\n")
} else {
  conn <- manc_neuprint()
  mn.ad <- manc_neuprint_meta("/type:MNad.*", conn = conn)
  en.ab <- manc_neuprint_meta("/type:EN.*",   conn = conn) |>
    filter(somaNeuromere %in% paste0("A", 1:8))
  targets <- bind_rows(
    mn.ad |> mutate(group = "MNad motor"),
    en.ab |> mutate(group = "abdominal EN")
  ) |>
    mutate(bodyid = as.character(bodyid)) |>
    select(bodyid, type, group, somaNeuromere, somaSide, exitNerve, predictedNt)
  saveRDS(targets, targets_f)
  cat("  fetched: nrow=", nrow(targets), "\n")
}
cat("  groups:\n"); print(table(targets$group))

cat("\nStep 2a: 3-D view of a sample of MANC abdominal neurons\n")
fig3d <- file.path(OUT_DIR, "abdominal_manc_neurons_3d.png")
if (!file.exists(fig3d)) {
  conn <- manc_neuprint()
  set.seed(1)
  sample_ids <- c(
    sample(targets$bodyid[targets$group == "MNad motor"],   8, replace = FALSE),
    sample(targets$bodyid[targets$group == "abdominal EN"], 6, replace = FALSE)
  )
  neurons <- manc_read_neurons(sample_ids, conn = conn)
  meta <- targets[match(names(neurons), targets$bodyid), ]
  cols <- ifelse(meta$group == "MNad motor", "#1f78b4", "#ff7f00")
  # Render off-screen via rgl + scene-graph.
  rgl::open3d()
  rgl::par3d(windowRect = c(50, 50, 1000, 800))
  if (requireNamespace("malevnc", quietly = TRUE)) {
    try(plot3d(MANC.surf, alpha = 0.05, col = "grey80"), silent = TRUE)
  }
  for (i in seq_along(neurons)) {
    plot3d(neurons[[i]], soma = 1500, lwd = 2, col = cols[i], add = TRUE)
  }
  rgl::view3d(theta = 0, phi = 0, fov = 0, zoom = 0.7)
  rgl::snapshot3d(fig3d, fmt = "png", webshot = FALSE, top = TRUE)
  rgl::close3d()
  cat("  wrote", fig3d, "\n")
} else {
  cat("  cached:", fig3d, "\n")
}

cat("\nStep 2b: NB v3_9_0 hits per body (cached)\n")
hits_f <- file.path(CACHE, "nb_hits_vnc.rds")
# Cap targets at 50 — neuronbridge_hits() uses plyr::rbind.fill
# internally, which is quadratic in result rows; on broad MANC hit
# lists (~thousands of rows per call) the full 273-body run hangs.
# 50 bodies is enough to populate a representative cell-type x KDRC
# heatmap for the figure.
TARGETS_CAP <- 50L
targets_run <- if (nrow(targets) > TARGETS_CAP) {
  set.seed(2)
  bind_rows(
    targets %>% filter(group == "MNad motor")  %>% slice_sample(n = floor(TARGETS_CAP * 0.7)),
    targets %>% filter(group == "abdominal EN") %>% slice_sample(n = ceiling(TARGETS_CAP * 0.3))
  )
} else {
  targets
}
cat("  capped to", nrow(targets_run), "bodies for the figure run\n")
if (file.exists(hits_f)) {
  hits.vnc <- readRDS(hits_f)
  cat("  cached: nrow=", nrow(hits.vnc), "\n")
} else {
  KEEP <- c("publishedName","libraryName","normalizedScore",
            "alignmentSpace","anatomicalArea","gender")
  mip.map <- list()
  for (bid in targets_run$bodyid) {
    info <- try(neuronbridge_info(bid, dataset = "by_body",
                                  version = NB_VERSION), silent = TRUE)
    if (inherits(info, "try-error") || is.null(info) || !nrow(info)) next
    manc <- info[grepl("MANC", info$libraryName, ignore.case = TRUE), ]
    if (!nrow(manc)) next
    manc$bodyid <- bid
    mip.map[[bid]] <- manc[1, c("bodyid", "nb.id", "libraryName")]
  }
  mip.map <- bind_rows(mip.map)
  cat("  resolved", nrow(mip.map), "/", nrow(targets_run), "MANC MIPs\n")

  hits <- list()
  pb <- txtProgressBar(min = 0, max = nrow(mip.map), style = 3)
  for (i in seq_len(nrow(mip.map))) {
    setTxtProgressBar(pb, i)
    h <- try(neuronbridge_hits(mip.map$nb.id[i], version = NB_VERSION),
             silent = TRUE)
    if (inherits(h, "try-error") || is.null(h) || !nrow(h)) next
    keep <- intersect(KEEP, colnames(h))
    h <- as.data.frame(h[, keep, drop = FALSE])
    h$normalizedScore <- as.numeric(h$normalizedScore)
    h$bodyid <- mip.map$bodyid[i]
    # best MIP per published line for this body
    h <- h[order(h$normalizedScore, decreasing = TRUE), , drop = FALSE]
    h <- h[!duplicated(h$publishedName), , drop = FALSE]
    hits[[i]] <- h
  }
  close(pb)
  hits <- bind_rows(hits)
  hits.vnc <- hits |>
    filter(grepl("FlyLight|Split", libraryName, ignore.case = TRUE),
           !grepl("Brain", libraryName, ignore.case = TRUE))
  hits.vnc <- merge(hits.vnc,
                    targets_run[, c("bodyid","type","group","somaNeuromere","exitNerve")],
                    by = "bodyid", all.x = TRUE)
  saveRDS(hits.vnc, hits_f)
  cat("  fetched: nrow=", nrow(hits.vnc), "\n")
}

cat("\nStep 2c: per-neuron elbow top-5 (drop frac 0.75, floor 0.20, cap 5)\n")
top_with_elbow <- function(df, cap = 5, drop_frac = 0.75, floor_frac = 0.20) {
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
hits.top <- hits.vnc |>
  group_by(bodyid) |>
  group_modify(~ top_with_elbow(.x)) |>
  ungroup()
candidate.lines <- sort(unique(hits.top$publishedName))
cat("  hits.top rows:", nrow(hits.top),
    "  unique lines:", length(candidate.lines), "\n")

cat("\nStep 2d: line-count summary plot (cell-type × line richness)\n")
# Per-cell-type summary: how many candidate lines per type, coloured by group.
ct_summary <- hits.top |>
  group_by(type, group) |>
  summarise(n_lines = n_distinct(publishedName),
            n_bodies = n_distinct(bodyid),
            best_score = max(normalizedScore, na.rm = TRUE),
            .groups = "drop") |>
  arrange(desc(n_lines))

p_lines <- ct_summary |>
  head(40) |>
  mutate(type = factor(type, levels = rev(type))) |>
  ggplot(aes(x = type, y = n_lines, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(`MNad motor` = "#1f78b4",
                               `abdominal EN` = "#ff7f00")) +
  labs(x = NULL, y = "n distinct VNC GAL4 lines (top-5 elbow)",
       title = "Candidate VNC GAL4 lines per MANC abdominal cell type",
       subtitle = "Top-40 cell types by line richness; from NB v3_9_0 colour-depth hits") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 7))
ggsave(file.path(OUT_DIR, "abdominal_lines_per_celltype.png"),
       p_lines, width = 8, height = 9, dpi = 300, bg = "white")
cat("  wrote abdominal_lines_per_celltype.png\n")

cat("\nStep 3: KDRC cross-reference\n")
if (SKIP_KDRC) {
  cat("  --skip-kdrc set — skipping live chromote runs\n")
  cat("  Synthesising a representative KDRC table from a sample of\n")
  cat("  candidate.lines so the heatmap demonstrates the layout the\n")
  cat("  live pipeline produces (no real biological inference!).\n")
  REGIONS <- c("Crop","PV","R1","R2","R3","R4","R5","MHJ","Hindgut","Rectum")
  # Sample 8 candidate lines and assign 1-3 random regions to each.
  set.seed(7)
  demo_lines <- sample(unique(hits.top$publishedName),
                       min(8, length(unique(hits.top$publishedName))))
  example_kdrc <- do.call(rbind, lapply(demo_lines, function(ln) {
    n <- sample.int(3, 1)
    tibble::tibble(publishedName = ln,
                   region = sample(REGIONS, n))
  }))
} else {
  ext_f <- file.path(CACHE, "kdrc_lookup.rds")
  if (file.exists(ext_f)) {
    example_kdrc <- readRDS(ext_f)
    cat("  KDRC cached: nrow=", nrow(example_kdrc), "\n")
  } else {
    cat("  Running chromote lookup (slow — ~1h per 500 lines)\n")
    REGIONS <- c("Crop","PV","R1","R2","R3","R4","R5","MHJ","Hindgut","Rectum")
    sess <- kdrc_start_session()
    kdrc <- kdrc_lookup_lines(candidate.lines, session = sess)
    kdrc_close_session(sess)
    halves <- split_halves(candidate.lines)
    extra.lines <- setdiff(unique(c(halves$ad, halves$dbd)),
                           c(candidate.lines, NA_character_))
    sess <- kdrc_start_session()
    kdrc_ext <- kdrc_lookup_lines(extra.lines, session = sess)
    kdrc_close_session(sess)
    to_long <- function(tab) {
      tab |>
        filter(kgut_hit == "Y") |>
        select(publishedName = query_line, all_of(REGIONS)) |>
        pivot_longer(all_of(REGIONS), names_to = "region", values_to = "flag") |>
        filter(flag == "Y") |>
        select(publishedName, region)
    }
    example_kdrc <- bind_rows(to_long(kdrc), to_long(kdrc_ext)) |> distinct()
    saveRDS(example_kdrc, ext_f)
    cat("  fetched: nrow=", nrow(example_kdrc), "\n")
  }
}

cat("\nStep 4: cell-type × KDRC region heatmap\n")
REGIONS <- c("Crop","PV","R1","R2","R3","R4","R5","MHJ","Hindgut","Rectum")
peripheral <- hits.top |>
  inner_join(example_kdrc, by = "publishedName") |>
  select(bodyid, type, group, publishedName, region)

if (!nrow(peripheral)) {
  warning("No peripheral hits at all — heatmap will be empty. Pass without --skip-kdrc to run live KDRC.")
}

mat.df <- peripheral |>
  count(type, region, name = "n") |>
  complete(type   = sort(unique(hits.top$type)),
           region = REGIONS, fill = list(n = 0))

# Filter to types with at least one hit so heatmap is informative.
keep_types <- mat.df |> group_by(type) |>
  summarise(s = sum(n), .groups = "drop") |> filter(s > 0) |> pull(type)
mat.df <- mat.df |> filter(type %in% keep_types)
gl <- hits.top |> distinct(type, group)
group_lookup <- setNames(gl$group, gl$type)

p_heat <- mat.df |>
  mutate(group = group_lookup[type],
         region = factor(region, levels = REGIONS)) |>
  ggplot(aes(x = region, y = reorder(type, n, sum), fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            colour = "white", size = 3) +
  scale_fill_gradient(low = "grey90", high = "#08306b", name = "n hits") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "KDRC gut region", y = NULL,
       title = "Abdominal MANC cell-type × KDRC gut-region matches",
       subtitle = ifelse(SKIP_KDRC,
                         "(demo data — pass without --skip-kdrc to run the live chromote lookup)",
                         "Direct hits + AD/DBD hemidriver hits combined")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 7))
ggsave(file.path(OUT_DIR, "abdominal_kdrc_heatmap.png"),
       p_heat, width = 7.5, height = 8, dpi = 300, bg = "white")
cat("  wrote abdominal_kdrc_heatmap.png\n")

cat("\nDone.\n")
