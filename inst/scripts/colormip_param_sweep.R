#!/usr/bin/env Rscript
# Tune `colormip_search()` parameters against an NBLAST ground-truth
# ranking. The default colormip parameters were tuned this way for the
# SREN query in the search_banc_with_lm vignette (Spearman 0.88 ->
# 0.901 vs NBLAST after disabling mirror and widening z_tolerance).
#
# Usage: edit the four paths under "Inputs" and run.
#
# Inputs
# ------
#   QUERY     a colour-depth MIP PNG (the query)
#   LIB       directory of candidate library MIP PNGs (one per neuron)
#   NB_SCORES an .rds with a data.frame: root_888 (chr) + nblast_mean
#             (numeric); the ground-truth ranking the sweep optimises
#             agreement with
#   META      optional .feather with cell metadata (root_888 + ...)
#             -- only used if you want a `cm_top6 ∩ nb_top6` overlap
#             count beyond the rank correlation
#
# Output
# ------
#   one .rds (`cmip_grid.rds`) with the full grid + Spearman + top-6
#   overlap, sorted by overlap then correlation. Best params for the
#   SREN query: threshold = 100, z_tolerance = 8, xy_shift = 3,
#   mirror = FALSE.

suppressMessages({
  library(neuronbridger)
})

# ----- Inputs (edit these) ----------------------------------------------
QUERY     <- "/path/to/query_mip.png"
LIB       <- "/path/to/library_dir"
NB_SCORES <- "/path/to/nblast_scores.rds"   # data.frame(root_888, nblast_mean)

# ----- Sweep grid -------------------------------------------------------
grid <- expand.grid(
  threshold   = c(60L, 100L, 200L, 500L),
  z_tolerance = c(0L, 1L, 2L, 4L, 6L, 8L),
  xy_shift    = c(0L, 2L, 3L, 4L),
  mirror      = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

# ----- Pool restricted to candidates NBLAST scored ----------------------
nb        <- readRDS(NB_SCORES)
all_paths <- list.files(LIB, "\\.png$", full.names = TRUE)
all_ids   <- sub(".*/(\\d+)_in_JRC2018.*", "\\1", all_paths)
LIB_FILES <- all_paths[all_ids %in% nb$root_888]
cat("Library size for grid:", length(LIB_FILES), "\n")
cat("Grid size:", nrow(grid), "\n")

run_one <- function(i) {
  g <- grid[i, ]
  res <- colormip_search(
    query       = QUERY,
    library     = LIB_FILES,
    threshold   = g$threshold,
    z_tolerance = g$z_tolerance,
    xy_shift    = g$xy_shift,
    mirror      = g$mirror,
    mc.cores    = max(1L, parallel::detectCores() - 1L),
    verbose     = FALSE)
  res$root_888 <- sub(".*/(\\d+)_in_JRC2018.*", "\\1", res$path)
  m <- merge(res[, c("root_888", "score")], nb,
             by = "root_888", sort = FALSE)
  rho <- suppressWarnings(
    stats::cor(m$score, m$nblast_mean, method = "spearman"))
  m_cm    <- m[order(-m$score), ]
  nb_top6 <- head(nb$root_888[order(-nb$nblast_mean)], 6L)
  cm_top6 <- head(m_cm$root_888, 6L)
  data.frame(
    threshold    = g$threshold,
    z_tolerance  = g$z_tolerance,
    xy_shift     = g$xy_shift,
    mirror       = g$mirror,
    spearman     = round(rho, 3),
    top6_overlap = length(intersect(cm_top6, nb_top6)))
}

out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  cat(sprintf("[%2d/%2d] thr=%d ztol=%d xy=%d mir=%s ... ",
              i, nrow(grid), grid$threshold[i], grid$z_tolerance[i],
              grid$xy_shift[i], grid$mirror[i]))
  r <- run_one(i)
  cat(sprintf("rho=%.3f  top6=%d\n", r$spearman, r$top6_overlap))
  r
}))

out <- out[order(-out$top6_overlap, -out$spearman), ]
saveRDS(out, "cmip_grid.rds")

cat("\n=== Top 10 (by overlap then rho) ===\n")
print(head(out, 10))
