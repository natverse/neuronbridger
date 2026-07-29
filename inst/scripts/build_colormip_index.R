#!/usr/bin/env Rscript
#
# build_colormip_index.R --- one-shot builder for the sparse library-index
# cache used by colormip_search(library_index = ...) /
# colormip_search_against_banc.R --cache <rds>.
#
# Reads every PNG in the BANC neuron colour-MIP library once, runs the
# depth-LUT indexer, and writes a single RDS holding per-candidate
# (fg_idx, z) pairs. Subsequent searches skip PNG decode + LUT indexing
# entirely, giving a ~5-10x per-candidate speedup on top of any
# query-side optimisations.
#
# Usage:
#   Rscript inst/scripts/build_colormip_index.R REGION [--threshold X]
#                                                        [--mc.cores N]
#                                                        [--data-root P]
#
#     REGION          "brain" | "VNC"
#     --threshold X   integer brightness cutoff to bake into the cache
#                     (must match the threshold passed to colormip_search
#                     later). Default 100L, same as colormip_search default.
#     --mc.cores N    parallel workers for the one-shot decode pass.
#                     Default 8L.
#     --data-root P   override the local data root. Default:
#                     ~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat
#                     (or $NEURONBRIDGER_DATA_ROOT).
#
# Local layout:
#   banc_colormips/{JRC2018_UNISEX_20x_HR,JRC2018_VNC_UNISEX_461}/
#     source library MIPs (populated once via `gsutil rsync`)
#   banc_colormips/<region>_cache_thr<T>.rds
#     output cache (large -- ~14 GB brain, ~3 GB VNC at typical sparsity)

suppressMessages({
  pkg_dir <- normalizePath(getwd(), mustWork = FALSE)
  if (file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      identical(read.dcf(file.path(pkg_dir, "DESCRIPTION"), "Package")[[1]],
                "neuronbridger")) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    library(neuronbridger)
  }
})

DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"
args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) && i < length(args)) return(args[i + 1L])
  default
}
strip_flag <- function(name) {
  i <- which(args == name)
  if (length(i)) args <<- args[-c(i, i + 1L)]
}

data_root <- get_flag("--data-root", Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT))
threshold <- as.integer(get_flag("--threshold", "100"))
mc_cores  <- as.integer(get_flag("--mc.cores", "8"))
strip_flag("--data-root"); strip_flag("--threshold"); strip_flag("--mc.cores")

if (!length(args)) stop("usage: build_colormip_index.R REGION [--threshold X] [--mc.cores N] [--data-root P]")
REGION <- match.arg(args[1], c("brain", "VNC"))
DATA_ROOT <- path.expand(data_root)
LIB_DIR <- file.path(DATA_ROOT, "banc_colormips",
                     if (REGION == "brain") "JRC2018_UNISEX_20x_HR"
                     else "JRC2018_VNC_UNISEX_461")
OUT_RDS <- file.path(DATA_ROOT, "banc_colormips",
                     sprintf("%s_cache_thr%d.rds", REGION, threshold))

if (!dir.exists(LIB_DIR))
  stop("library dir not found: ", LIB_DIR,
       "\nRun: gsutil -m rsync -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_colormips/template_alignment_240721/",
       if (REGION == "brain") "JRC2018_UNISEX_20x_HR/" else "JRC2018_VNC_UNISEX_461/",
       " ", LIB_DIR, "/")

n <- length(list.files(LIB_DIR, "\\.png$"))
cat(sprintf("== building %s cache: %d candidates, threshold=%d, mc.cores=%d ==\n",
            REGION, n, threshold, mc_cores))
cat(sprintf("  from: %s\n  to:   %s\n\n", LIB_DIR, OUT_RDS))

build_colormip_index(LIB_DIR, OUT_RDS,
                     threshold = threshold, mc.cores = mc_cores, verbose = TRUE)
