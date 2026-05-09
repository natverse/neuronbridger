#!/usr/bin/env Rscript
# Combine per-sample <name>_hits.csv files in ~/deng_to_banc/colormip_hits/
# into a single master CSV. Reads banc_id as character to preserve the
# 18-char IDs.

suppressMessages({
  library(readr); library(dplyr); library(tibble)
})

HITS_DIR <- path.expand("~/deng_to_banc/colormip_hits")
OUT_CSV  <- path.expand("~/deng_to_banc/colormip_hits/master.csv")

files <- list.files(HITS_DIR, pattern = "_hits\\.csv$", full.names = TRUE)
files <- setdiff(files, OUT_CSV)
cat("found", length(files), "per-sample CSVs\n")

all <- lapply(files, function(f) {
  readr::read_csv(f, show_col_types = FALSE,
                  col_types = readr::cols(banc_id = "c"))
})
master <- dplyr::bind_rows(all)
master <- dplyr::arrange(master, deng_sample, dplyr::desc(cm_score))
readr::write_csv(master, OUT_CSV)
cat(sprintf("wrote %d rows (%d samples) -> %s\n",
            nrow(master),
            length(unique(master$deng_sample)),
            OUT_CSV))
print(master |> dplyr::group_by(deng_sample) |>
              dplyr::summarise(n_hits = dplyr::n(),
                               n_nblast = sum(!is.na(nblast)),
                               max_cm = max(cm_score),
                               max_nb = if (any(!is.na(nblast))) max(nblast, na.rm = TRUE) else NA_real_))
