#!/usr/bin/env Rscript
#
# make_colormip_library.R --- generate NeuronBridge-compatible colour-MIPs
# from Deng/Kondo LM source volumes.
#
# Why this is separate from the *_to_banc.R registration scripts: those
# scripts target BANC voxel space (for visual EM overlay). NeuronBridge
# searches need MIPs on the canonical template grids:
#   brain : JRC2018_UNISEX_20x_HR  (1210 x 566 x 174 @ 0.519/0.519/1.0 um)
#   VNC   : JRC2018_VNC_UNISEX_HR  (573 x 1119 x 219 @ 0.461/0.461/0.7 um)
# `nrrd_to_mip(target_space = "VNC")` adds the standard 90-px black header
# that every Janelia VNC colour-MIP carries, so the output of this script
# is directly searchable against the BANC neuron MIP library at
#   gs://lee-lab.../neuron_colormips/template_alignment_240721/
#
# Per-sample work:
#   Stage 1   re-warp LSM (Deng) or NRRD (Kondo) -> NB target template
#             via Elastix (Deng) or the existing brain image-mode warp
#             chain (Kondo).
#   Stage 2   nrrd_to_mip(method = "direct") -> PNG.
#
# Usage:
#   Rscript inst/scripts/make_colormip_library.R DATASET REGION \
#     [test|all] [--force] [--data-root <path>]
#
#     DATASET         "deng" | "kondo"
#     REGION          "brain" | "VNC"  (Kondo is brain-only)
#     test            4-sample probe;  all = full cohort
#     --force         recompute even if the PNG already exists
#                     (default skips existing outputs).
#     --data-root P   override local data root. Default:
#                     ~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat
#                     (or $NEURONBRIDGER_DATA_ROOT).
#     --threshold X   MIP threshold passed to nrrd_to_mip(). Default 0.999
#                     (top 0.1% of nonzero voxels). Matches BANC library
#                     MIP sparsity so colormip_search actually finds real
#                     matches; the old "auto" (Triangle) default left LM
#                     stains ~40-80% foreground vs the ~0.001% typical of
#                     library MIPs, drowning the search signal.
#
# Local layout expected under <data-root>:
#   imaging-CCT-Bowen/         Deng raw LSMs
#   kondo_et_al_2020/nrrd/     Kondo IS2 NRRDs
#   templates/                 JRC2018_UNISEX_20x_HR.nrrd, JRC2018_VNC_UNISEX_HR.nrrd
#   {deng,kondo}_to_banc/      outputs (this script writes
#                              colormip_library/{brain,VNC}/<name>.png here)

suppressMessages({
  pkg_dir <- normalizePath(getwd(), mustWork = FALSE)
  if (file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      identical(read.dcf(file.path(pkg_dir, "DESCRIPTION"), "Package")[[1]],
                "neuronbridger")) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    library(neuronbridger)
  }
  library(nat)
})

# --- knobs ---------------------------------------------------------
DEFAULT_DATA_ROOT <- "~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"

args <- commandArgs(trailingOnly = TRUE)
dr_idx <- which(args == "--data-root")
data_root <- if (length(dr_idx) && dr_idx < length(args)) args[dr_idx + 1L] else Sys.getenv("NEURONBRIDGER_DATA_ROOT", unset = DEFAULT_DATA_ROOT)
DATA_ROOT <- path.expand(data_root)
if (length(dr_idx)) args <- args[-c(dr_idx, dr_idx + 1L)]

force <- "--force" %in% args
args  <- args[args != "--force"]

thr_idx <- which(args == "--threshold")
mip_threshold <- if (length(thr_idx) && thr_idx < length(args)) {
  as.numeric(args[thr_idx + 1L])
} else 0.999
if (length(thr_idx)) args <- args[-c(thr_idx, thr_idx + 1L)]

if (length(args) < 2) {
  stop("usage: make_colormip_library.R DATASET REGION [test|all] [--force] [--threshold X] [--data-root <path>]")
}
dataset <- match.arg(args[1], c("deng", "kondo"))
region  <- match.arg(args[2], c("brain", "VNC"))
mode    <- if (length(args) >= 3) args[3] else "test"
if (dataset == "kondo" && region == "VNC")
  stop("Kondo et al. 2020 is brain-only.")

DENG_LSM_DIR <- file.path(DATA_ROOT, "imaging-CCT-Bowen")
KONDO_DIR    <- file.path(DATA_ROOT, "kondo_et_al_2020", "nrrd")
TPL <- list(
  brain = file.path(DATA_ROOT, "templates", "JRC2018_UNISEX_20x_HR.nrrd"),
  VNC   = file.path(DATA_ROOT, "templates", "JRC2018_VNC_UNISEX_HR.nrrd")
)
OUT_ROOT <- file.path(DATA_ROOT,
                      sprintf("%s_to_banc", dataset),
                      "colormip_library", region)
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
cat("DATA_ROOT: ", DATA_ROOT, "\n", sep = "")

# Disable the brain-specific central-brain mask — it's on the JRC2018U grid,
# not JRC2018U_HR or JRC2018VNCU_HR.
options(neuronbridger.jrc2018u_fixed_mask = "")

# --- Deng LSM source ----------------------------------------------
parse_lsm <- function(fname) {
  bn <- sub("\\.lsm$", "", basename(fname), ignore.case = TRUE)
  m <- regmatches(bn, regexec("^(.*?)[-_].*?(BRAIN|VNC)([-_][FM])?", bn,
                              ignore.case = TRUE))[[1]]
  if (!length(m) || length(m) < 3L) return(NULL)
  list(gene = m[[2]], region = toupper(m[[3]]),
       sex  = if (length(m) >= 4L && nzchar(m[[4]]))
                toupper(sub("[-_]", "", m[[4]])) else "U")
}

list_deng <- function(region, mode) {
  pat <- if (region == "brain") "BRAIN.*\\.lsm$" else "VNC.*\\.lsm$"
  src <- list.files(DENG_LSM_DIR, pattern = pat, full.names = TRUE,
                    recursive = FALSE, ignore.case = TRUE)
  if (mode == "test") {
    head(src[grepl("CapaR", basename(src), ignore.case = TRUE)], 4)
  } else src
}

list_kondo <- function(mode) {
  src <- list.files(KONDO_DIR, pattern = "_02_warp.*\\.nrrd$",
                    full.names = TRUE)
  if (mode == "test") head(src, 3) else src
}

# --- per-sample worker --------------------------------------------
process_deng <- function(lsm) {
  meta <- parse_lsm(lsm); if (is.null(meta)) return(NULL)
  if ((region == "brain" && meta$region != "BRAIN") ||
      (region == "VNC"   && meta$region != "VNC"))   return(NULL)
  name <- sprintf("%s_%s-%s", meta$gene, meta$region, meta$sex)
  out_png <- file.path(OUT_ROOT, paste0(name, ".png"))
  if (file.exists(out_png) && !force) {
    message("SKIP ", name, " (PNG exists; pass --force to recompute)")
    return(out_png)
  }
  message("\n--- ", name, " ---")
  raw_dir <- file.path(OUT_ROOT, "raw")
  reg_dir <- file.path(OUT_ROOT, "elastix", name)
  chans <- lsm_to_nrrd(lsm, raw_dir, prefix = name, verbose = FALSE)
  reg <- lm_to_jrc2018u_elastix(
    nc82 = chans[["nc82"]], gfp = chans[["gfp"]],
    output_dir = reg_dir, template = TPL[[region]],
    threads = 8L, verbose = FALSE
  )
  v <- nat::read.nrrd(reg$gfp_jrc2018u); storage.mode(v) <- "integer"
  png_path <- nrrd_to_mip(v, save = TRUE, format = "png",
                          target_space = region, method = "direct",
                          savefolder = file.path(OUT_ROOT, "_mip_tmp"),
                          threshold = mip_threshold)
  file.rename(png_path, out_png)
  unlink(file.path(OUT_ROOT, "_mip_tmp"), recursive = TRUE)
  # Keep the warped GFP NRRD next to the PNG; the colormip search step
  # needs it for LM dotprops + voxel-attribution.
  warp_dst <- file.path(OUT_ROOT, paste0(name, ".nrrd"))
  if (file.exists(reg$gfp_jrc2018u)) file.copy(reg$gfp_jrc2018u, warp_dst,
                                                overwrite = TRUE)
  unlink(reg_dir, recursive = TRUE)
  message("  wrote ", out_png)
  out_png
}

process_kondo <- function(nrrd) {
  bn <- sub("\\.nrrd$", "", basename(nrrd), ignore.case = TRUE)
  m  <- regmatches(bn, regexec("^IS2_(.+?)_no(\\d+)_", bn))[[1]]
  if (length(m) < 3) { message("can't parse ", bn); return(NULL) }
  name <- sprintf("%s_no%s", m[[2]], m[[3]])
  out_png <- file.path(OUT_ROOT, paste0(name, ".png"))
  if (file.exists(out_png)) { message("SKIP ", name); return(out_png) }
  message("\n--- ", name, " ---")
  # Use the existing IS2->FCWB->JRC2018F image-mode chain, then resample
  # to JRC2018U_HR. This is the same Stage A+B path lsm_to_banc_layer
  # uses minus the BANC step.
  fcwb_nrrd <- file.path(OUT_ROOT, paste0(name, "_FCWB.nrrd"))
  jrcf_nrrd <- file.path(OUT_ROOT, paste0(name, "_JRC2018F.nrrd"))
  jrcu_nrrd <- file.path(OUT_ROOT, paste0(name, "_JRC2018U_HR.nrrd"))
  # IS2 -> FCWB (CMTK)
  is2_to_fcwb_cmtk(nrrd, fcwb_nrrd, verbose = FALSE)
  # FCWB -> JRC2018F (H5)
  fcwb_to_jrc2018f_h5(fcwb_nrrd, jrcf_nrrd, threads = 8L, verbose = FALSE)
  # JRC2018F -> JRC2018U_HR via nat::resample to the HR grid.
  v_jrcf <- nat::read.nrrd(jrcf_nrrd)
  hr <- nat::read.nrrd(TPL[[region]])
  v_hr <- nat::resample(v_jrcf, hr)   # nat dispatches to im3d resample
  nat::write.nrrd(v_hr, jrcu_nrrd, dtype = "byte")
  v <- nat::read.nrrd(jrcu_nrrd); storage.mode(v) <- "integer"
  png_path <- nrrd_to_mip(v, save = TRUE, format = "png",
                          target_space = region, method = "direct",
                          savefolder = file.path(OUT_ROOT, "_mip_tmp"),
                          threshold = mip_threshold)
  file.rename(png_path, out_png)
  unlink(file.path(OUT_ROOT, "_mip_tmp"), recursive = TRUE)
  # Keep <name>.nrrd next to the PNG for downstream search NBLAST.
  warp_dst <- file.path(OUT_ROOT, paste0(name, ".nrrd"))
  if (file.exists(jrcu_nrrd)) file.rename(jrcu_nrrd, warp_dst)
  for (f in c(fcwb_nrrd, jrcf_nrrd)) if (file.exists(f)) file.remove(f)
  message("  wrote ", out_png)
  out_png
}

# --- main ----------------------------------------------------------
if (dataset == "deng") {
  targets <- list_deng(region, mode)
  cat(sprintf("== %d Deng %s LSMs (mode '%s') ==\n",
              length(targets), region, mode))
  out <- lapply(targets, function(s) tryCatch(process_deng(s),
                                              error = function(e) {
                                                message("FAILED: ", conditionMessage(e)); NULL
                                              }))
} else {
  targets <- list_kondo(mode)
  cat(sprintf("== %d Kondo brain NRRDs (mode '%s') ==\n",
              length(targets), mode))
  out <- lapply(targets, function(s) tryCatch(process_kondo(s),
                                              error = function(e) {
                                                message("FAILED: ", conditionMessage(e)); NULL
                                              }))
}
n_ok <- sum(!vapply(out, is.null, logical(1)))
cat(sprintf("\nwrote %d / %d PNGs to %s\n", n_ok, length(targets), OUT_ROOT))
