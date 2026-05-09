#!/usr/bin/env Rscript
#
# deng_colormip_search.R --- Per-sample colour-MIP + NBLAST against BANC.
#
# For every Deng et al. 2019 LM stain in the local mirror:
#   Stage 1  re-warp the LSM to JRC2018U_HR (brain) or
#            JRC2018VNCU_HR (VNC). Output: <name>_in_<template>.nrrd.
#   Stage 2  render a colour-depth MIP (PNG) via `nrrd_to_mip()`.
#   Stage 3  search the BANC neuron MIP library at
#            gs://lee-lab_brain-and-nerve-cord-fly-connectome/
#                                  neuron_colormips/template_alignment_240721/
#            for top-K hits via `colormip_search()`.
#   Stage 4  NBLAST top-K BANC neurons against the LM-derived dotprops.
#   Stage 5  emit one CSV per sample + one master CSV at the end.
#
# Resumable: skips samples whose master CSV row already exists.
#
# Usage:
#   Rscript inst/scripts/deng_colormip_search.R [test|brain|vnc|all] [TOP_K]
#     test  -- 4 CapaR samples (BRAIN-{F,M} + VNC-{F,M})
#     brain -- every BRAIN file
#     vnc   -- every VNC file
#     all   -- both
#     TOP_K -- per-sample top hits to keep (default 25)

suppressMessages({
  # Prefer the working-tree neuronbridger if present at NB_DEV_ROOT or
  # the cwd's parent path (the script ships under inst/scripts/).
  dev_root <- Sys.getenv("NB_DEV_ROOT", unset = "")
  if (!nzchar(dev_root)) {
    cwd <- normalizePath(getwd(), mustWork = FALSE)
    if (file.exists(file.path(cwd, "DESCRIPTION")) &&
        identical(read.dcf(file.path(cwd, "DESCRIPTION"), "Package")[[1]],
                  "neuronbridger")) dev_root <- cwd
  }
  if (nzchar(dev_root) && dir.exists(dev_root)) {
    devtools::load_all(dev_root, quiet = TRUE)
  } else {
    library(neuronbridger)
  }
  library(bancr)
  library(nat)
  library(nat.flybrains)
  library(nat.templatebrains)
  library(nat.jrcbrains)
  library(nat.nblast)
  nat.jrcbrains::register_saalfeldlab_registrations()
  library(tibble); library(dplyr); library(readr)
})

LSM_DIR  <- "/Users/asbates/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/imaging-CCT-Bowen"
RAW_DIR  <- "~/deng_to_banc/raw"
WORK_DIR <- "~/deng_to_banc/colormip_search"
LIB_DIR  <- "~/banc_colormips"
OUT_CSV  <- "~/deng_to_banc/colormip_hits/master.csv"

WORK_DIR <- path.expand(WORK_DIR)
LIB_DIR  <- path.expand(LIB_DIR)
OUT_CSV  <- path.expand(OUT_CSV)
dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1) args[1] else "test"
TOP_K <- if (length(args) >= 2) as.integer(args[2]) else 25L

# --- locate samples (mirrors deng_to_banc.R parser) -------------------
parse_lsm <- function(fname) {
  bn <- sub("\\.lsm$", "", basename(fname), ignore.case = TRUE)
  m  <- regmatches(bn, regexec("^(.*?)[-_].*?(BRAIN|VNC)([-_][FM])?", bn,
                               ignore.case = TRUE))[[1]]
  if (!length(m) || length(m) < 3L) return(NULL)
  list(gene   = m[[2]],
       region = toupper(m[[3]]),
       sex    = if (length(m) >= 4L && nzchar(m[[4]])) toupper(sub("[-_]", "", m[[4]])) else "U")
}

GENE_TEST <- c("CapaR-L-T1-BRAIN-F", "CapaR-L-T1-BRAIN-M",
               "CapaR-L-T1-VNC-F",   "CapaR-L-T1-VNC-M")

list_lsms <- function(mode) {
  if (mode == "test") {
    sapply(GENE_TEST, function(s) file.path(LSM_DIR, paste0(s, ".lsm")))
  } else {
    pat <- switch(mode,
                  brain = "BRAIN.*\\.lsm$",
                  vnc   = "VNC.*\\.lsm$",
                  all   = "(BRAIN|VNC).*\\.lsm$",
                  stop("mode must be test|brain|vnc|all"))
    list.files(LSM_DIR, pattern = pat, full.names = TRUE,
               recursive = FALSE, ignore.case = TRUE)
  }
}

samples <- list_lsms(mode)
samples <- samples[file.exists(samples)]
cat(sprintf("found %d LSMs in mode '%s'\n", length(samples), mode))

# --- per-sample pipeline ---------------------------------------------

process_one <- function(lsm, top_k = TOP_K) {
  meta <- parse_lsm(lsm)
  if (is.null(meta)) {
    message("  parse failed: ", basename(lsm)); return(NULL)
  }
  region <- if (meta$region == "VNC") "VNC" else "brain"
  bn     <- sub("\\.lsm$", "", basename(lsm), ignore.case = TRUE)
  name   <- sprintf("%s_%s-%s",
                    sub("[-_].*", "", bn, ignore.case = TRUE),
                    meta$region, meta$sex)
  out_dir   <- file.path(WORK_DIR, name)
  out_csv   <- file.path(dirname(OUT_CSV), paste0(name, "_hits.csv"))
  if (file.exists(out_csv)) {
    message("  SKIP (have CSV): ", name)
    return(readr::read_csv(out_csv, show_col_types = FALSE,
                           col_types = readr::cols(banc_id = "c")))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Stage 0: extract NC82 + GFP NRRDs (re-use prior raw if present)
  raw_dir <- path.expand(RAW_DIR)
  nc82 <- file.path(raw_dir, paste0(bn, "_nc82.nrrd"))
  gfp  <- file.path(raw_dir, paste0(bn, "_gfp.nrrd"))
  if (!file.exists(nc82) || !file.exists(gfp)) {
    chans <- lsm_to_nrrd(lsm, raw_dir, prefix = bn, verbose = FALSE)
    nc82  <- chans[["nc82"]]; gfp <- chans[["gfp"]]
  }

  # Stage 1: re-warp to NeuronBridge target template
  template <- switch(region,
                     brain = "~/templates/JRC2018_UNISEX_20x_HR.nrrd",
                     VNC   = "~/templates/JRC2018_VNC_UNISEX_HR.nrrd")
  options(neuronbridger.jrc2018u_fixed_mask = "")
  reg <- lm_to_jrc2018u_elastix(
    nc82       = nc82,
    gfp        = gfp,
    output_dir = file.path(out_dir, "elastix"),
    template   = path.expand(template),
    threads    = 8L,
    verbose    = FALSE
  )

  # Stage 2: colour-MIP from warped GFP volume
  v <- nat::read.nrrd(reg$gfp_jrc2018u)
  storage.mode(v) <- "integer"
  mip_path <- nrrd_to_mip(
    v,
    save = TRUE, format = "png",
    target_space = region, method = "direct",
    savefolder = file.path(out_dir, "color_mips")
  )

  # Stage 3: search BANC library
  lib_subdir <- if (region == "brain") "JRC2018_UNISEX_20x_HR" else "JRC2018_VNC_UNISEX_461"
  lib <- file.path(LIB_DIR, lib_subdir)
  if (!dir.exists(lib))
    stop("BANC ", region, " MIP library not synced: ", lib,
         "\n  rsync from gs://lee-lab.../neuron_colormips/template_alignment_240721/", lib_subdir, "/")
  hits <- colormip_search(
    query     = mip_path,
    library   = lib,
    threshold = 100L, z_tolerance = 2L, xy_shift = 2L,
    mirror    = TRUE, top_n = top_k,
    mc.cores  = 8L, verbose = FALSE
  )
  if (!nrow(hits)) return(NULL)

  # Stage 4: NBLAST top-K vs LM dotprops
  hit_ids <- sub(".*/(\\d+)_in_.*", "\\1", hits$path)
  fg <- which(v > 30, arr.ind = TRUE)
  if (nrow(fg) > 50000L) fg <- fg[sample.int(nrow(fg), 50000L), , drop = FALSE]
  vd <- diag(attr(v, "header")[["space directions"]])
  lm_pts <- sweep(fg - 1L, 2, vd, "*")
  lm_dp  <- nat::dotprops(lm_pts, k = 5L)

  # Per-neuron NBLAST — skip ones that fail to load or warp.
  target_tb <- if (region == "brain") "JRC2018U"  else "JRCVNC2018U"
  source_tb <- if (region == "brain") "JRC2018F"  else "JRCVNC2018F"
  region_tps <- if (region == "brain") "brain" else "vnc"
  nblast_scores <- vapply(seq_along(hit_ids), function(i) {
    tryCatch({
      m <- bancr::banc_read_neuron_meshes(hit_ids[i])
      if (length(m) == 0L) return(NA_real_)
      m_jf <- bancr::banc_to_JRC2018F(m,
                                      method = "tpsreg",
                                      banc.units = "nm",
                                      region = region_tps)
      m_tg <- nat.templatebrains::xform_brain(m_jf,
                                              sample = source_tb, reference = target_tb)
      pts <- nat::xyzmatrix(m_tg)
      if (nrow(pts) < 6L) return(NA_real_)
      dp  <- nat::dotprops(pts, k = 5L)
      as.numeric(nat.nblast::nblast(dp, target = lm_dp, normalised = TRUE))
    }, error = function(e) NA_real_)
  }, numeric(1))
  cat(sprintf("  NBLAST: %d / %d scored\n",
              sum(!is.na(nblast_scores)), length(nblast_scores)))

  out <- tibble::tibble(
    deng_sample = name,
    region      = region,
    banc_id     = as.character(hit_ids),
    cm_score    = hits$score,
    cm_n_match  = hits$n_match,
    cm_dx       = hits$dx,
    cm_dy       = hits$dy,
    cm_mirror   = hits$mirror,
    nblast      = nblast_scores
  )
  readr::write_csv(out, out_csv)
  out
}

# --- main loop --------------------------------------------------------

all_hits <- list()
for (i in seq_along(samples)) {
  cat(sprintf("\n--- [%d/%d] %s ---\n", i, length(samples), basename(samples[i])))
  t0 <- Sys.time()
  out <- tryCatch(process_one(samples[i]),
                  error = function(e) {
                    message("  FAILED: ", conditionMessage(e)); NULL
                  })
  if (!is.null(out)) all_hits[[length(all_hits) + 1]] <- out
  cat(sprintf("  done in %.0fs\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if (length(all_hits)) {
  master <- dplyr::bind_rows(all_hits)
  readr::write_csv(master, OUT_CSV)
  cat(sprintf("\nwrote %d rows to %s\n", nrow(master), OUT_CSV))
} else {
  cat("\nno samples produced hits\n")
}
