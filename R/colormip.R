#' @title Make a NeuronBridge-style colour-depth MIP from an .nrrd image volume
#'
#' @description Encode the most-anterior occupied z-slice at each \code{(x, y)}
#' pixel as a colour drawn from a 256-entry depth lookup table (blue =
#' anterior, red = posterior; PsychedelicRainBow2-like ramp), and write the
#' resulting RGB MIP to a \code{.png} or \code{.tif}. Three back-ends are
#' offered, all targeting the same Janelia ColorMIP / Color Depth MIP
#' specification (Otsuna \emph{et al.} 2018,
#' \href{https://www.biorxiv.org/content/10.1101/318006}{bioRxiv}):
#'
#' \describe{
#'   \item{\code{method = "direct"} (default)}{Pure R; vectorised
#'     \code{which.max} along z, LUT lookup, transpose. No external
#'     dependencies beyond \code{nat}, \code{png} and (for TIFF) \code{tiff}.
#'     Pixel-perfect against the BANC Python port.}
#'   \item{\code{method = "python"}}{Calls the Python implementation by
#'     Stephan Gerhard distributed in
#'     \href{https://github.com/jasper-tms/the-BANC-fly-connectome/blob/main/fanc/render_neurons.py}{the
#'     BANC connectome package} (\code{fanc.render_neurons.make_colormip}) via
#'     \code{reticulate}. The R helper loads each NRRD, hands the volume to
#'     numpy, and runs the inner colour-MIP algorithm with the BANC
#'     \code{depth_lut}. Use this to validate the R port byte-for-byte.}
#'   \item{\code{method = "fiji"}}{Launches FIJI/ImageJ and runs Janelia's
#'     \code{Color_Depth_MIP_batch_0308_2021.ijm} macro (the canonical
#'     reference implementation). Requires a working FIJI install with the
#'     ColorMIP plugins (\href{https://github.com/JaneliaSciComp/ColorMIP_Mask_Search}{search tool},
#'     \href{https://github.com/JaneliaSciComp/ColorMIP_Mask_Search/tree/master/ColorDepthMIP_Generator}{MIP generator}).
#'     The macro asks for input and output folders interactively; the
#'     \code{input}/\code{savefolder} arguments are not used in this path.}
#' }
#'
#' @param input either a path to a single \code{.nrrd} file, a path to a
#' folder of \code{.nrrd} files (each will be processed), or an \code{im3d}/3D
#' numeric array. Ignored when \code{method = "fiji"}.
#' @param savefolder folder in which to write the MIP image(s). Defaults to a
#' \code{color_mips/} folder next to the input. Required when \code{input} is
#' an in-memory array and \code{save = TRUE}. Ignored when \code{method = "fiji"}.
#' @param method which back-end to use. One of \code{"direct"} (default,
#' pure R), \code{"python"} (BANC Python via \code{reticulate}), or
#' \code{"fiji"} (Janelia FIJI macro).
#' @param target_space which template space the input volume is in.
#' \code{"brain"} (default) corresponds to **\code{JRC2018U_HR}** (a.k.a.
#' \code{JRC2018_UNISEX_20x_HR}, dims \code{c(1210, 566, 174)}, voxdims
#' \code{c(0.5189, 0.5189, 1.0)}) — the NeuronBridge brain space.
#' \code{"VNC"} corresponds to **\code{JRC2018VNCU_HR}** (a.k.a.
#' \code{JRC2018_VNC_UNISEX_461}, dims \code{c(573, 1119, 219)}, voxdims
#' \code{c(0.461, 0.461, 0.7)}) — the NeuronBridge VNC space. Both are
#' the high-resolution variants Janelia's pre-computed ColorMIPs are
#' served at. For \code{"VNC"} a 90-pixel black header is prepended so
#' output is dimensionally identical to Janelia VNC colormips.
#' @param threshold one of \code{"auto"} (default), \code{"triangle"},
#' \code{"otsu"}, \code{"none"}, or a numeric value. Controls how
#' foreground voxels are separated from background:
#' \itemize{
#'   \item \code{"auto"} keeps the legacy \code{> 0} mask if the input
#'     looks binary (≤ 2 distinct values, e.g. a synthetic neuron
#'     skeleton from \code{\link{root_id_to_nrrd}}); otherwise applies
#'     Triangle thresholding to suppress background staining in
#'     real-world LM data.
#'   \item \code{"triangle"} (Zack \emph{et al.} 1977): draws a line
#'     from the histogram peak to the highest-value bin and picks the
#'     bin of maximum perpendicular distance. Robust on histograms with
#'     a dominant background peak and a long signal tail (typical
#'     confocal data).
#'   \item \code{"otsu"} (1979): maximises between-class variance,
#'     ignoring zero voxels by default. Works on bimodal foregrounds.
#'   \item \code{"none"} keeps every voxel \code{> 0} (the original
#'     behaviour; correct for binary inputs).
#'   \item A \strong{numeric in [0, 1]} is interpreted as a quantile
#'     threshold against the non-zero voxels (\code{0.99} = top 1\%).
#'   \item A \strong{numeric > 1} is a raw intensity cutoff.
#' }
#' @param denoise one of \code{"auto"} (default), \code{"median3d"} or
#' \code{"none"}. Controls a pre-MIP denoising step that runs before
#' \code{threshold}: \code{"median3d"} applies a 3 × 3 × 3 median
#' filter to suppress salt-and-pepper noise (requires the
#' \pkg{mmand} package). \code{"auto"} enables \code{"median3d"} for
#' grayscale input and skips it for binary input. The 3-D filter is
#' equivalent to FIJI's \code{Mask Median Subtraction} idiom that the
#' Janelia ColorMIP macro uses for the same purpose.
#' @param format output image format. One of \code{"png"} (default) or \code{"tiff"}.
#' @param save logical; if \code{TRUE} (default) write the MIP to disk and
#' return the file path(s) invisibly. If \code{FALSE}, return the MIP as a
#' \code{height x width x 3} numeric array with values in \code{[0, 1]}.
#' @param overwrite logical; if \code{FALSE} (default) skip inputs whose MIP
#' file already exists.
#' @param fiji.path \emph{(method = "fiji")} path to FIJI executable. If
#' \code{NULL}, looks up \code{getOption("jimpipeline.fiji")}, then \code{$PATH},
#' then \code{/Applications/Fiji.app/...}.
#' @param macro \emph{(method = "fiji")} path to the FIJI macro. Defaults to
#' \code{Color_Depth_MIP_batch_0308_2021.ijm} (shipped with this package and
#' expected to be installed under FIJI's \code{plugins/Macros}).
#' @param MinMem,MaxMem \emph{(method = "fiji")} JVM memory limits.
#' @param python \emph{(method = "python")} optional path to a Python
#' executable or virtualenv. If \code{NULL}, \code{reticulate}'s default
#' Python is used. If the BANC \code{depth_lut} cannot be imported from
#' \code{fanc.render_neurons}, the function falls back to an internal
#' verbatim copy of the same LUT.
#'
#' @return Invisibly, the path(s) to the written MIP file(s); or, when
#' \code{save = FALSE}, the MIP as a \code{height x width x 3} array (or list
#' of such, when processing a folder). For \code{method = "fiji"}, returns
#' the result of the macro call.
#'
#' @examples
#' \dontrun{
#' # Render a flywire neuron into a JRC2018U_HR-sized .nrrd, then color-MIP it
#' # in pure R (no FIJI, no Python):
#' root_id_to_nrrd("720575940632295751",
#'                 reference  = "JRC2018U_HR",
#'                 savefolder = "~/asta_sez_mips")
#' nrrd_to_mip("~/asta_sez_mips", target_space = "brain")
#'
#' # Validate against Stephan Gerhard's Python port:
#' nrrd_to_mip("~/asta_sez_mips", method = "python")
#'
#' # Or hand off to the FIJI macro (interactive folder picker):
#' nrrd_to_mip(method = "fiji",
#'             fiji.path = "/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx")
#' }
#' @seealso \code{\link{root_id_to_nrrd}}, \code{\link{ngl_scene_to_nrrd}}
#' @export
nrrd_to_mip <- function(input = NULL,
                        savefolder = NULL,
                        method = c("direct", "python", "fiji"),
                        target_space = c("brain", "VNC"),
                        threshold = "auto",
                        denoise = "auto",
                        format = c("png", "tiff"),
                        save = TRUE,
                        overwrite = FALSE,
                        fiji.path = NULL,
                        macro = NULL,
                        MinMem = "2500m",
                        MaxMem = "2500m",
                        python = NULL) {
  method <- match.arg(method)
  target_space <- match.arg(target_space)
  format <- match.arg(format)

  if (method == "fiji") {
    if (!is.null(input))
      message("method = 'fiji' opens the macro's interactive folder picker; ",
              "`input` and `savefolder` are not forwarded.")
    return(nrrd_to_mip_fiji_impl(fiji.path = fiji.path, macro = macro,
                                 MinMem = MinMem, MaxMem = MaxMem))
  }

  fn <- if (method == "direct") {
    function(vol) colormip_from_array(vol, target_space = target_space,
                                      threshold = threshold, denoise = denoise)
  } else {
    function(vol) colormip_from_array_python(vol, target_space = target_space,
                                             threshold = threshold,
                                             denoise = denoise,
                                             python = python)
  }

  process_one <- function(vol, outfile) {
    cmip <- fn(vol)
    if (!save) return(cmip)
    if (format == "png") {
      png::writePNG(cmip, target = outfile)
    } else {
      if (!requireNamespace("tiff", quietly = TRUE))
        stop("Install the 'tiff' package to write TIFF output, ",
             "or use format = 'png'.")
      tiff::writeTIFF(cmip, where = outfile, bits.per.sample = 8L)
    }
    message("Wrote: ", outfile)
    invisible(outfile)
  }

  # Folder dispatch
  if (is.character(input) && length(input) == 1 && dir.exists(input)) {
    files <- list.files(input, pattern = "\\.nrrd$", full.names = TRUE,
                        ignore.case = TRUE)
    if (!length(files)) stop("No .nrrd files found in: ", input)
    if (is.null(savefolder)) savefolder <- file.path(input, "color_mips")
    dir.create(savefolder, showWarnings = FALSE, recursive = TRUE)
    out <- lapply(files, function(f) {
      nrrd_to_mip(f, savefolder = savefolder, method = method,
                  target_space = target_space,
                  threshold = threshold, denoise = denoise,
                  format = format,
                  save = save, overwrite = overwrite, python = python)
    })
    names(out) <- tools::file_path_sans_ext(basename(files))
    return(invisible(out))
  }

  # Single-file dispatch
  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    if (!grepl("\\.nrrd$", input, ignore.case = TRUE))
      stop("Expected an .nrrd file: ", input)
    if (is.null(savefolder)) savefolder <- file.path(dirname(input), "color_mips")
    base <- tools::file_path_sans_ext(basename(input))
    outfile <- file.path(savefolder, paste0(base, "_colormip.", format))
    if (save) {
      dir.create(savefolder, showWarnings = FALSE, recursive = TRUE)
      if (!overwrite && file.exists(outfile)) {
        message("Skipping (exists): ", outfile)
        return(invisible(outfile))
      }
    }
    vol <- nat::read.nrrd(input)
    return(process_one(vol, outfile))
  }

  # In-memory array dispatch
  if (is.array(input) || is.numeric(input)) {
    if (length(dim(input)) != 3)
      stop("Expected a 3D volume; got dim = ", paste(dim(input), collapse = "x"))
    if (save && is.null(savefolder))
      stop("`savefolder` must be supplied when `input` is an in-memory array, ",
           "or set `save = FALSE` to return the MIP.")
    outfile <- if (save) {
      dir.create(savefolder, showWarnings = FALSE, recursive = TRUE)
      file.path(savefolder, paste0("colormip.", format))
    } else NULL
    return(process_one(input, outfile))
  }

  stop("`input` must be a path to an .nrrd file, a folder of .nrrd files, ",
       "or a 3D array (got: ", class(input)[1], ").")
}

# Core algorithm in pure R: 3D array (x, y, z) -> 2D MIP (height, width, 3) in
# [0, 1]. Mirrors fanc.render_neurons.make_colormip (Stephan Gerhard, BANC),
# with optional pre-MIP denoising + thresholding for noisy LM volumes.
colormip_from_array <- function(vol,
                                target_space = c("brain", "VNC"),
                                threshold    = "auto",
                                denoise      = "auto") {
  target_space <- match.arg(target_space)
  vol <- colormip_preprocess(vol, threshold = threshold, denoise = denoise)
  d <- dim(vol)
  nx <- d[1]; ny <- d[2]; nz <- d[3]

  # Most-anterior occupied slice per (x, y); first occurrence on a binary
  # input matches np.argmax in the BANC port.
  zmax <- apply(vol, c(1, 2), which.max)            # 1..nz
  maxv <- apply(vol, c(1, 2), max)                  # 0 where no signal
  fg <- maxv > 0

  z255 <- as.integer((zmax - 1) / (nz - 1) * 255)   # 0..255, truncated
  z255[!fg] <- 0L

  lut <- colormip_depth_lut                         # 256 x 3 integer matrix
  rgb_flat <- matrix(0L, nrow = nx * ny, ncol = 3)
  fg_flat  <- as.vector(fg)
  idx_flat <- as.vector(z255) + 1L                  # 1..256
  rgb_flat[fg_flat, ] <- lut[idx_flat[fg_flat], ]

  cmip <- array(rgb_flat, dim = c(nx, ny, 3))
  cmip <- aperm(cmip, c(2, 1, 3))                   # rows = y, cols = x

  if (target_space == "VNC") {
    header <- array(0L, dim = c(90L, dim(cmip)[2], 3L))
    cmip <- abind3(header, cmip)
  }

  cmip / 255
}

# Apply optional denoise + threshold to a 3-D volume before colour-MIP.
# Returns the same shape as the input. Always integer-valued.
colormip_preprocess <- function(vol, threshold = "auto", denoise = "auto") {
  # Cheap binarity probe: sample up to 5000 voxels (skip any zero); the
  # data is binary iff the sampled non-zero voxels are all equal. The
  # synthetic neuron renderings produced by root_id_to_nrrd() / as.im3d
  # are typically {0, 255}; LM volumes have hundreds of distinct values.
  if (length(vol) > 5000L) {
    s <- vol[seq.int(1L, length(vol), length.out = 5000L)]
  } else s <- vol
  s <- s[s > 0]
  is_binary <- !length(s) || length(unique.default(s)) == 1L

  # --- Denoise ---
  use_denoise <- if (identical(denoise, "auto")) {
    if (is_binary) "none" else "median3d"
  } else if (is.character(denoise) && length(denoise) == 1L) {
    match.arg(denoise, c("none", "median3d"))
  } else stop("`denoise` must be 'auto', 'none' or 'median3d'.")
  if (use_denoise == "median3d") {
    if (!requireNamespace("mmand", quietly = TRUE))
      stop("`denoise = \"median3d\"` requires the 'mmand' package: ",
           "install.packages('mmand').")
    k <- mmand::shapeKernel(c(3L, 3L, 3L), type = "box")
    vol <- mmand::medianFilter(vol, k)
  }

  # --- Threshold ---
  cutoff <- if (is.numeric(threshold) && length(threshold) == 1L &&
                !is.na(threshold)) {
    nzv <- vol[vol > 0]
    if (threshold > 0 && threshold < 1 && length(nzv))
      as.numeric(stats::quantile(nzv, threshold, names = FALSE))
    else as.numeric(threshold)
  } else if (is.character(threshold) && length(threshold) == 1L) {
    method <- if (identical(threshold, "auto")) {
      if (is_binary) "none" else "triangle"
    } else match.arg(threshold, c("auto", "triangle", "otsu", "none"))
    switch(method,
      none     = 0,
      triangle = colormip_triangle_threshold(vol),
      otsu     = colormip_otsu_threshold(vol))
  } else stop("`threshold` must be 'auto', 'triangle', 'otsu', 'none', ",
              "or a single numeric.")

  if (cutoff > 0) vol[vol <= cutoff] <- 0L
  storage.mode(vol) <- "integer"
  vol
}

# Triangle thresholding (Zack et al. 1977). Operates on the histogram of
# integer voxel values 0..nbins-1. Robust on histograms with a dominant
# background peak and a long signal tail.
colormip_triangle_threshold <- function(vol, nbins = 256L) {
  h <- tabulate(as.integer(vol) + 1L, nbins = nbins)
  pk <- which.max(h)
  last <- max(which(h > 0))
  if (last <= pk + 1L) return(pk - 1L)
  xs <- as.numeric(pk:last); ys <- as.numeric(h[pk:last])
  x1 <- xs[1]; y1 <- ys[1]; x2 <- xs[length(xs)]; y2 <- ys[length(ys)]
  dist <- abs((y2 - y1) * xs - (x2 - x1) * ys + x2 * y1 - y2 * x1)
  as.integer(xs[which.max(dist)] - 1L)
}

# Otsu thresholding (1979). By default ignores zero voxels so the
# threshold separates background staining from real signal rather than
# zero-padding from everything else.
colormip_otsu_threshold <- function(vol, nbins = 256L, ignore_zero = TRUE) {
  vv <- if (ignore_zero) vol[vol > 0] else as.integer(vol)
  if (!length(vv)) return(0L)
  h <- tabulate(as.integer(vv) + 1L, nbins = nbins) / length(vv)
  cs <- cumsum(h); m <- cumsum((seq_len(nbins) - 1L) * h)
  mt <- m[nbins]
  num <- (mt * cs - m)^2
  den <- cs * (1 - cs); den[den == 0] <- NA
  as.integer(which.max(num / den) - 1L)
}

# Python back-end: hand the volume to numpy and run the colormip math from
# fanc.render_neurons.make_colormip. Falls back to an internal copy of
# depth_lut if the BANC package is not installed.
colormip_from_array_python <- function(vol,
                                       target_space = c("brain", "VNC"),
                                       threshold    = "auto",
                                       denoise      = "auto",
                                       python = NULL) {
  target_space <- match.arg(target_space)
  vol <- colormip_preprocess(vol, threshold = threshold, denoise = denoise)
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("method = 'python' requires the 'reticulate' package. Install with ",
         "install.packages('reticulate').")
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)

  # The PyPI `banc` package exposes `banc.render_neurons`; the source repo
  # ships the same module at `fanc.render_neurons`. Try both before falling
  # back to the bundled LUT.
  banc_lut <- NULL
  for (mod in c("banc.render_neurons", "fanc.render_neurons")) {
    banc_lut <- tryCatch(reticulate::import(mod, convert = TRUE)$depth_lut,
                         error = function(e) NULL)
    if (!is.null(banc_lut)) break
  }
  if (is.null(banc_lut)) {
    message("BANC python package not importable; using bundled depth_lut.")
    banc_lut <- lapply(seq_len(nrow(colormip_depth_lut)),
                       function(i) as.integer(colormip_depth_lut[i, ]))
  }

  np <- reticulate::import("numpy", convert = TRUE)
  py <- reticulate::py
  py[["_colormip_vol"]] <- np$asarray(reticulate::r_to_py(vol), dtype = "uint8")
  py[["_colormip_lut"]] <- banc_lut
  py[["_colormip_vnc"]] <- target_space == "VNC"

  reticulate::py_run_string("
import numpy as np
from skimage.color import rgb2hsv, hsv2rgb
_vol = _colormip_vol
_lut = list(_colormip_lut)
target_zmax = _vol.shape[-1] - 1
image_zmaxes = np.argmax(_vol, axis=2).astype(np.float64)
image_zmaxes = (image_zmaxes / target_zmax * 255).astype(np.uint32)
colormip = np.zeros((_vol.shape[0], _vol.shape[1], 3), dtype=np.uint8)
for i in range(image_zmaxes.shape[0]):
    for j in range(image_zmaxes.shape[1]):
        if image_zmaxes[i, j] != 0:
            hsv = rgb2hsv(np.array(_lut[int(image_zmaxes[i, j])], dtype=np.double))
            hsv[2] = 255
            colormip[i, j, :] = hsv2rgb(hsv).astype(np.uint8)
colormip = colormip.transpose(1, 0, 2).astype(np.uint8)
if _colormip_vnc:
    colormip = np.vstack([np.zeros((90, colormip.shape[1], 3), dtype=np.uint8),
                          colormip]).astype(np.uint8)
")

  reticulate::py$colormip / 255
}

# FIJI back-end: launches FIJI and runs Janelia's macro. Replaces the
# pre-refactor standalone nrrd_to_mip().
nrrd_to_mip_fiji_impl <- function(fiji.path = NULL,
                                  macro = NULL,
                                  MinMem = "2500m",
                                  MaxMem = "2500m") {
  fiji.path <- fiji(fijiPath = fiji.path)
  message("fiji: ", fiji.path)
  if (is.null(macro)) {
    os <- get_os()
    if (os == "windows") {
      macro <- normalizePath(file.path(
        dirname(fiji.path),
        "plugins\\Macros\\Color_Depth_MIP_batch_0308_2021.ijm"))
    } else {
      macro <- normalizePath(
        "/Applications/Fiji.app/plugins/Macros/Color_Depth_MIP_batch_0308_2021.ijm")
    }
  } else {
    macro <- normalizePath(macro)
  }
  runFijiMacro(
    macro = macro, macroArg = "", headless = FALSE, batch = FALSE,
    MinMem = MinMem, MaxMem = MaxMem, IncrementalGC = TRUE,
    Threads = NULL, fijiArgs = NULL, javaArgs = NULL, ijArgs = NULL,
    fijiPath = fiji.path, DryRun = FALSE
  )
}

# Tiny rbind for 3-D arrays along the first axis (no abind dependency)
abind3 <- function(a, b) {
  da <- dim(a); db <- dim(b)
  if (!identical(da[2:3], db[2:3]))
    stop("abind3: dims 2 and 3 must match")
  out <- array(0, dim = c(da[1] + db[1], da[2], da[3]))
  out[seq_len(da[1]), , ] <- a
  out[da[1] + seq_len(db[1]), , ] <- b
  out
}

# 256-entry depth LUT used by the Janelia ColorMIP algorithm
# (PsychedelicRainBow2-like ramp: blue = anterior, red = posterior).
# Verbatim from fanc.render_neurons.depth_lut.
colormip_depth_lut <- matrix(c(
  127,0,255, 125,3,255, 124,6,255, 122,9,255, 121,12,255, 120,15,255,
  119,18,255, 118,21,255, 116,24,255, 115,27,255, 114,30,255, 113,33,255,
  112,36,255, 110,39,255, 109,42,255, 108,45,255, 106,48,255, 105,51,255,
  104,54,255, 103,57,255, 101,60,255, 100,63,255, 99,66,255, 98,69,255,
  96,72,255, 95,75,255, 94,78,255, 93,81,255, 92,84,255, 90,87,255, 89,90,255,
  87,93,255, 86,96,255, 84,99,255, 83,102,255, 81,105,255, 80,108,255,
  78,111,255, 77,114,255, 75,117,255, 74,120,255, 72,123,255, 71,126,255,
  69,129,255, 68,132,255, 66,135,255, 65,138,255, 63,141,255, 62,144,255,
  60,147,255, 59,150,255, 57,153,255, 56,156,255, 54,159,255, 53,162,255,
  51,165,255, 50,168,255, 48,171,255, 47,174,255, 45,177,255, 44,180,255,
  42,183,255, 41,186,255, 39,189,255, 38,192,255, 36,195,255, 35,198,255,
  33,201,255, 32,204,255, 30,207,255, 29,210,255, 27,213,255, 26,216,255,
  24,219,255, 23,222,255, 21,225,255, 20,228,255, 18,231,255, 16,234,255,
  14,237,255, 12,240,255, 9,243,255, 6,246,255, 3,249,255, 1,252,255,
  0,254,255, 3,255,252, 6,255,249, 9,255,246, 12,255,243, 15,255,240,
  18,255,237, 21,255,234, 24,255,231, 27,255,228, 30,255,225, 33,255,222,
  36,255,219, 39,255,216, 42,255,213, 45,255,210, 48,255,207, 51,255,204,
  54,255,201, 57,255,198, 60,255,195, 63,255,192, 66,255,189, 69,255,186,
  72,255,183, 75,255,180, 78,255,177, 81,255,174, 84,255,171, 87,255,168,
  90,255,165, 93,255,162, 96,255,159, 99,255,156, 102,255,153, 105,255,150,
  108,255,147, 111,255,144, 114,255,141, 117,255,138, 120,255,135, 123,255,132,
  126,255,129, 129,255,126, 132,255,123, 135,255,120, 138,255,117, 141,255,114,
  144,255,111, 147,255,108, 150,255,105, 153,255,102, 156,255,99, 159,255,96,
  162,255,93, 165,255,90, 168,255,87, 171,255,84, 174,255,81, 177,255,78,
  180,255,75, 183,255,72, 186,255,69, 189,255,66, 192,255,63, 195,255,60,
  198,255,57, 201,255,54, 204,255,51, 207,255,48, 210,255,45, 213,255,42,
  216,255,39, 219,255,36, 222,255,33, 225,255,30, 228,255,27, 231,255,24,
  234,255,21, 237,255,18, 240,255,15, 243,255,12, 246,255,9, 249,255,6,
  252,255,3, 254,255,0, 255,252,3, 255,249,6, 255,246,9, 255,243,12,
  255,240,15, 255,237,18, 255,234,21, 255,231,24, 255,228,27, 255,225,30,
  255,222,33, 255,219,36, 255,216,39, 255,213,42, 255,210,45, 255,207,48,
  255,204,51, 255,201,54, 255,198,57, 255,195,60, 255,192,63, 255,189,66,
  255,186,69, 255,183,72, 255,180,75, 255,177,78, 255,174,81, 255,171,84,
  255,168,87, 255,165,90, 255,162,93, 255,159,96, 255,156,99, 255,153,102,
  255,150,105, 255,147,108, 255,144,111, 255,141,114, 255,138,117, 255,135,120,
  255,132,123, 255,129,126, 255,126,129, 255,123,132, 255,120,135, 255,117,138,
  255,114,141, 255,111,144, 255,108,147, 255,105,150, 255,102,153, 255,99,156,
  255,96,159, 255,93,162, 255,90,165, 255,87,168, 255,84,171, 255,81,173,
  255,78,174, 255,75,175, 255,72,176, 255,69,177, 255,66,178, 255,63,179,
  255,60,180, 255,57,181, 255,54,182, 255,51,183, 255,48,184, 255,45,185,
  255,42,186, 255,39,187, 255,36,188, 255,33,189, 255,30,190, 255,27,191,
  255,24,192, 255,21,193, 255,18,194, 255,15,195, 255,12,196, 255,9,197,
  255,6,198, 255,3,199, 255,0,200
), ncol = 3, byrow = TRUE)
