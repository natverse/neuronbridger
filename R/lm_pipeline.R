#' @title Bridge a source-space LM volume into BANC voxel coordinates
#'
#' @description Three-stage \strong{image-mode} pipeline that lands a
#' light-microscopy volume on the BANC EM voxel grid:
#'
#' \enumerate{
#'   \item \strong{IS2 -> FCWB} via CMTK
#'   \code{reformatx} (image-mode resampling using the bridging
#'   registration shipped in \pkg{nat.flybrains}).
#'   \item \strong{FCWB -> JRC2018F} via the saalfeldlab
#'   \code{template-building} \code{RenderTransformed} Java tool
#'   driven by the H5 displacement field
#'   (\code{nat.jrcbrains/JRC2018F_FCWB.h5}).
#'   \item \strong{JRC2018F -> BANC voxel space} via Elastix
#'   \code{transformix} using \pkg{bancr}'s bundled
#'   \code{brain_240721/BANC_to_template.txt} transform chain.
#' }
#'
#' Output lands on the BANC precomputed grid
#' (\code{2400 x 924 x 789} at 400 nm) ready for
#' \code{\link{nrrd_to_precomputed}}.
#'
#' Unlike the earlier points-mode prototype, every stage runs in
#' native code on full image volumes. Per-volume runtime on a recent
#' Apple-silicon Mac with 8 threads:
#'
#' \itemize{
#'   \item Stage A (CMTK reformatx): ~30 s
#'   \item Stage B (RenderTransformed): ~3 min
#'   \item Stage C (transformix):     ~1 min
#'   \item Stage D (nrrd_to_precomputed): ~30 s
#'   \item \strong{total}: ~5 min/volume
#' }
#'
#' @section External dependencies:
#' \itemize{
#'   \item CMTK \code{reformatx} on \code{PATH}, or set
#'         \code{options(neuronbridger.cmtk_bin = "/path/to/cmtk-bin")}.
#'   \item Elastix 5.x \code{transformix} on \code{PATH}.
#'   \item Java 17+ runtime.
#'   \item The saalfeldlab \code{RenderTransformed} fat JAR. By
#'         default looked up at
#'         \code{~/.local/share/neuronbridger/render-transformed.jar};
#'         override with the env var \code{NEURONBRIDGER_XFORMIMAGE_JAR}
#'         or \code{options(neuronbridger.xformimage_jar = "/abs/path/to.jar")}.
#'   \item \pkg{nat.jrcbrains} H5 registrations downloaded via
#'         \code{nat.jrcbrains::download_saalfeldlab_registrations()}.
#' }
#'
#' To build the JAR locally:
#' \preformatted{
#' git clone https://github.com/saalfeldlab/template-building.git
#' cd template-building
#' # remove the unresolvable org.itc dependency (RenderTransformed
#' # does not use it)
#' sed -i '' '/<groupId>org.itc<\/groupId>/,/<\/dependency>/d' pom.xml
#' mvn package -DskipTests
#' mkdir -p ~/.local/share/neuronbridger
#' cp target/application-jar-with-dependencies.jar \\
#'    ~/.local/share/neuronbridger/render-transformed.jar
#' }
#'
#' @section Status:
#' Hidden / internal helper while we benchmark across the Kondo
#' glutamate-receptor library. If stable, propose upstream as a thin
#' \code{xformimage.h5reg} wrapper for \pkg{nat.h5reg}.
#'
#' @name lm_pipeline
#' @keywords internal
NULL


#' @rdname lm_pipeline
#'
#' @param input source-space NRRD path (e.g. an IS2-aligned Kondo
#'   stack).
#' @param output target output NRRD path on the FCWB grid.
#' @param sample_space source template name. Currently only
#'   \code{"IS2"} is wired up (the only Kondo source space we ship a
#'   target-grid + bridging registration for here). Patches welcome.
#' @param right_shift integer; right-shift applied to the source
#'   volume before warping (default \code{4} = 12-bit -> 8-bit, which
#'   matches Kondo et al.'s 12-bit-packed-in-uint16 stacks). The
#'   downstream Elastix output is uint8, so without this conversion
#'   12-bit signal saturates at 255 across the brain. Pass \code{0}
#'   for data already in the 0-255 range.
#' @param verbose logical; emit stage timing.
#' @return path to the FCWB-aligned NRRD.
#' @export
is2_to_fcwb_cmtk <- function(input,
                             output,
                             sample_space = "IS2",
                             right_shift  = 4L,
                             verbose      = TRUE) {
  stopifnot(file.exists(input), nzchar(output))
  if (!identical(sample_space, "IS2"))
    stop("Only sample_space = 'IS2' is currently supported.")

  reformatx <- .neuronbridger_reformatx()
  cmtk_reg  <- system.file("extdata", "bridgingregistrations",
                           "FCWB_IS2.list",
                           package = "nat.flybrains")
  if (!nzchar(cmtk_reg) || !dir.exists(cmtk_reg))
    stop("nat.flybrains FCWB_IS2.list not found.")

  # 12-bit -> 8-bit conversion before warping. Kondo NRRDs are 16-bit
  # files holding 12-bit data; the Elastix Stage C clips at 255, so
  # forwarding 16-bit through the chain saturates everything bright.
  src_for_cmtk <- input
  tmp_8bit     <- NULL
  if (is.numeric(right_shift) && right_shift > 0) {
    v <- nat::read.nrrd(input)
    hdr <- attr(v, "header")
    voxdims_um <- diag(hdr[["space directions"]])
    cap <- (2L ^ (8L + as.integer(right_shift))) - 1L
    v8 <- as.integer(pmin(pmax(as.integer(v), 0L), cap) %/%
                       (2L ^ as.integer(right_shift)))
    dim(v8) <- dim(v)
    v_im <- nat::im3d(v8, dims = dim(v),
                      voxdims = voxdims_um,
                      origin  = c(0, 0, 0))
    tmp_8bit <- tempfile(pattern = "is2_8bit_", fileext = ".nrrd")
    nat::write.nrrd(v_im, tmp_8bit, dtype = "byte", enc = "gzip")
    src_for_cmtk <- tmp_8bit
  }
  on.exit(if (!is.null(tmp_8bit) && file.exists(tmp_8bit)) file.remove(tmp_8bit),
          add = TRUE)

  # FCWB grid: 1769 x 1026 x 108 at 0.318967, 0.318427, 1.0 um.
  # CMTK's parser balks on >7 decimal-place voxel sizes, so 6 dp.
  target_grid <- "1769,1026,108:0.318967,0.318427,1.0"

  t0 <- Sys.time()
  rval <- system2(reformatx,
                  args = c("--target-grid", shQuote(target_grid),
                           "--floating",    shQuote(src_for_cmtk),
                           "--outfile",     shQuote(output),
                           shQuote(cmtk_reg)),
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("CMTK reformatx exited ", rval, " on ", basename(input))
  if (!file.exists(output))
    stop("reformatx claimed success but did not write ", output)

  if (verbose)
    message(sprintf("  [stage A] %s (%.1fs)",
                    basename(output),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(output)
}


#' @rdname lm_pipeline
#'
#' @param input FCWB-aligned NRRD path (output of
#'   \code{is2_to_fcwb_cmtk}).
#' @param output target output NRRD path on the JRC2018F grid.
#' @param h5_path path to a saalfeldlab H5 displacement field
#'   bridging FCWB -> JRC2018F. Defaults to
#'   \code{nat.jrcbrains}'s \code{JRC2018F_FCWB.h5}.
#' @param threads integer; number of threads passed to the JAR
#'   (\code{-q}). Default \code{8}.
#' @return path to the JRC2018F-aligned NRRD.
#' @export
fcwb_to_jrc2018f_h5 <- function(input,
                                output,
                                h5_path = NULL,
                                threads = 8L,
                                verbose = TRUE) {
  stopifnot(file.exists(input), nzchar(output))
  jar  <- .neuronbridger_xformimage_jar()
  java <- .neuronbridger_java_bin()
  if (is.null(h5_path)) {
    candidates <- c(
      "~/Library/Application Support/R/nat.jrcbrains/JRC2018F_FCWB/JRC2018F_FCWB.h5",
      "~/.local/share/R/nat.jrcbrains/JRC2018F_FCWB/JRC2018F_FCWB.h5"
    )
    h5_path <- Find(file.exists, sapply(candidates, path.expand))
    if (is.null(h5_path))
      stop("JRC2018F_FCWB.h5 not found. Run ",
           "nat.jrcbrains::download_saalfeldlab_registrations().")
  }

  # JRC2018F output grid: 1652 x 768 x 479 voxels at 0.38 um isotropic.
  # Saalfeld H5 dfield in JRC2018F_FCWB.h5 maps FCWB -> JRC2018F (forward).
  # For image resampling we need OUTPUT (JRC2018F) -> SOURCE (FCWB), i.e. the
  # inverse of dfield. RenderTransformed's `-i` flag inverts the next transform,
  # so we pass `-i <h5>:0/dfield:0/invdfield` to get JRC2018F -> FCWB.
  interval <- "0,0,0:1651,767,478"
  res      <- "0.38,0.38,0.38"
  xfm_arg  <- paste0(h5_path, ":0/dfield:0/invdfield")

  t0 <- Sys.time()
  rval <- system2(java,
                  args = c("-Xmx16g",
                           "-cp", shQuote(jar),
                           "process.RenderTransformed",
                           shQuote(input),
                           shQuote(output),
                           shQuote(interval),
                           "-r", res,
                           "-q", as.integer(threads),
                           "-i", shQuote(xfm_arg)),
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("RenderTransformed exited ", rval, " on ", basename(input))
  if (!file.exists(output))
    stop("RenderTransformed claimed success but did not write ", output)

  if (verbose)
    message(sprintf("  [stage B] %s (%.1fs)",
                    basename(output),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(output)
}


#' @rdname lm_pipeline
#'
#' @param input JRC2018F-aligned NRRD path.
#' @param output target output NRRD path on the BANC grid (uint8
#'   gzipped NRRD).
#' @param transform_dir directory holding the BANC Elastix transform
#'   chain (\code{BANC_to_template.txt} + \code{0/1/2/3_*.txt}).
#'   Defaults to \pkg{bancr}'s bundled \code{brain_240721/}.
#' @param transformix path to the \code{transformix} binary
#'   (default \code{"transformix"}, expected on \code{PATH}).
#' @return path to the BANC-aligned NRRD.
#' @export
jrc2018f_to_banc_elastix <- function(input,
                                     output,
                                     transform_dir = NULL,
                                     transformix   = "transformix",
                                     verbose       = TRUE) {
  input  <- path.expand(input)
  output <- path.expand(output)
  stopifnot(file.exists(input), nzchar(output))
  if (is.null(transform_dir)) {
    if (!requireNamespace("bancr", quietly = TRUE))
      stop("Either install bancr (which bundles brain_240721/) or pass `transform_dir`.")
    transform_dir <- system.file("extdata", "brain_240721", package = "bancr")
  }
  tp <- file.path(transform_dir, "BANC_to_template.txt")
  if (!file.exists(tp))
    stop("BANC_to_template.txt not found in ", transform_dir)
  if (Sys.which(transformix) == "")
    stop("`", transformix, "` not on PATH. Install Elastix 5.x.")

  t0 <- Sys.time()
  out_dir <- tempfile("transformix_")
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rval <- system2(transformix,
                  args = c("-in",  shQuote(input),
                           "-out", shQuote(out_dir),
                           "-tp",  shQuote(tp)),
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("transformix exited ", rval, " — see ", out_dir, "/transformix.log")
  raw <- file.path(out_dir, "result.nrrd")
  if (!file.exists(raw))
    stop("transformix did not write result.nrrd in ", out_dir)

  v <- nat::read.nrrd(raw)
  v[v < 0.5] <- 0
  v[v > 255] <- 255
  v <- as.integer(round(v))
  dim(v) <- c(2400L, 924L, 789L)
  v_im <- nat::im3d(v, dims = c(2400L, 924L, 789L),
                    voxdims = c(0.4, 0.4, 0.4),
                    origin  = c(0, 0, 0))
  nat::write.nrrd(v_im, output, dtype = "byte", enc = "gzip")

  if (verbose)
    message(sprintf("  [stage C] %s (%.1fs)",
                    basename(output),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(output)
}


#' @rdname lm_pipeline
#'
#' @param input source-space NRRD path.
#' @param gene short gene / receptor / peptide name (e.g.
#'   \code{"GluRIIA"}). Forms the precomputed dir name as
#'   \code{{gene}_{sample}}.
#' @param sample sample identifier (e.g. \code{"no1"} for Kondo
#'   first-specimen samples; \code{"no1"} / \code{"no2"} are the two
#'   prepared brains in Kondo et al. 2020).
#' @param channel source-stack channel string (default \code{"02"} =
#'   GFP / receptor channel for Kondo data; \code{"01"} would be the
#'   NC82 anatomical reference). Recorded in the registry only; the
#'   warping pipeline itself ignores this string.
#' @param dataset master folder under \code{light_level/} (default
#'   \code{"kondo_et_al_2020"}).
#' @param output_dir parent directory for intermediates + the final
#'   \code{{gene}_{sample}.ng/} precomputed dir.
#' @param source_path optional original storage path on the upstream
#'   repository; recorded in the registry entry.
#' @param chunk_size length-3 integer chunk size for the precomputed
#'   layer. Default \code{c(64, 64, 64)} matches the public BANC
#'   atlas.
#' @param keep_intermediates logical; keep the FCWB / JRC2018F / BANC
#'   intermediate NRRDs for inspection (default \code{FALSE}).
#' @param sample_space source template name (default \code{"IS2"}).
#' @param right_shift integer right-shift to apply to the source
#'   volume in Stage A (default \code{4} = 12-bit -> 8-bit).
#' @param threads passed to \code{fcwb_to_jrc2018f_h5}.
#' @return list with \code{precomputed_dir} (filesystem path),
#'   \code{registry_entry} (one-row tibble), and \code{timings}
#'   (named numeric, seconds per stage).
#' @export
lm_to_banc_layer <- function(input,
                             gene,
                             sample,
                             channel             = "02",
                             dataset             = "kondo_et_al_2020",
                             output_dir,
                             source_path         = NULL,
                             chunk_size          = c(64L, 64L, 64L),
                             keep_intermediates  = FALSE,
                             sample_space        = "IS2",
                             right_shift         = 4L,
                             threads             = 8L,
                             verbose             = TRUE) {
  stopifnot(file.exists(input),
            nzchar(gene), nzchar(sample), nzchar(dataset))
  if (missing(output_dir) || is.null(output_dir))
    stop("`output_dir` is required.")
  output_dir <- path.expand(output_dir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  name      <- paste(gene, sample, sep = "_")
  fcwb_nrrd <- file.path(output_dir, paste0(name, "_in_FCWB.nrrd"))
  jrcf_nrrd <- file.path(output_dir, paste0(name, "_in_JRC2018F.nrrd"))
  banc_nrrd <- file.path(output_dir, paste0(name, "_aligned240721_to_BANC.nrrd"))
  pc_dir    <- file.path(output_dir, paste0(name, ".ng"))

  t_all <- Sys.time()
  tA <- system.time(is2_to_fcwb_cmtk(input, fcwb_nrrd,
                                     sample_space = sample_space,
                                     right_shift  = right_shift,
                                     verbose = verbose))[["elapsed"]]
  tB <- system.time(fcwb_to_jrc2018f_h5(fcwb_nrrd, jrcf_nrrd,
                                        threads = threads,
                                        verbose = verbose))[["elapsed"]]
  tC <- system.time(jrc2018f_to_banc_elastix(jrcf_nrrd, banc_nrrd,
                                             verbose = verbose))[["elapsed"]]
  tD <- system.time(nrrd_to_precomputed(
    input      = banc_nrrd,
    output     = pc_dir,
    resolution = c(400, 400, 400),
    data_type  = "uint8",
    encoding   = "raw",
    chunk_size = chunk_size,
    overwrite  = TRUE
  ))[["elapsed"]]

  if (!keep_intermediates) {
    for (f in c(fcwb_nrrd, jrcf_nrrd, banc_nrrd))
      if (file.exists(f)) file.remove(f)
  }

  registry_entry <- tibble::tibble(
    dataset         = dataset,
    gene            = gene,
    sample          = sample,
    channel         = channel,
    name            = name,
    gs_url          = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/%s/%s.ng/",
                              dataset, name),
    alignment_space = "BANC",
    voxdims_nm      = list(c(400L, 400L, 400L)),
    source_filename = basename(input),
    source_path     = source_path %||% input,
    uploaded        = format(Sys.Date())
  )

  list(precomputed_dir = pc_dir,
       registry_entry  = registry_entry,
       timings         = c(stageA_is2_to_fcwb     = tA,
                           stageB_fcwb_to_jrcf    = tB,
                           stageC_jrcf_to_banc    = tC,
                           stageD_precomputed     = tD,
                           total                  = as.numeric(
                             difftime(Sys.time(), t_all, units = "secs"))))
}


# --- locator helpers ---------------------------------------------------

.neuronbridger_xformimage_jar <- function() {
  jar <- Sys.getenv("NEURONBRIDGER_XFORMIMAGE_JAR", unset = "")
  if (!nzchar(jar))
    jar <- getOption("neuronbridger.xformimage_jar",
                     default = path.expand("~/.local/share/neuronbridger/render-transformed.jar"))
  if (!file.exists(jar))
    stop("RenderTransformed JAR not found at: ", jar, "\n",
         "Build from saalfeldlab/template-building (see ?lm_pipeline) ",
         "or set NEURONBRIDGER_XFORMIMAGE_JAR / ",
         "options(neuronbridger.xformimage_jar = ...).")
  jar
}

.neuronbridger_java_bin <- function() {
  java <- Sys.getenv("NEURONBRIDGER_JAVA", unset = "")
  if (!nzchar(java)) java <- getOption("neuronbridger.java", default = "java")
  if (Sys.which(java) == "" && !file.exists(java))
    stop("Java runtime not found (tried: ", java,
         "). Install JDK 17+ and ensure `java` is on PATH or set ",
         "NEURONBRIDGER_JAVA / options(neuronbridger.java = ...).")
  java
}

.neuronbridger_reformatx <- function() {
  bin <- getOption("neuronbridger.cmtk_bin", default = "")
  exe <- if (nzchar(bin)) file.path(bin, "reformatx") else "reformatx"
  if (Sys.which(exe) == "" && !file.exists(exe)) {
    # try the homebrew cmtk-contrib location
    cand <- Sys.glob("/opt/homebrew/Cellar/cmtk*/bin/reformatx")
    if (length(cand)) return(cand[[1]])
    stop("CMTK `reformatx` not found. Install cmtk and ensure it's on PATH ",
         "or set options(neuronbridger.cmtk_bin = '/path/to/cmtk-bin').")
  }
  exe
}

`%||%` <- function(a, b) if (is.null(a) || (is.character(a) && !nzchar(a))) b else a


# ----------------------------------------------------------------------
# Deng et al. 2019 LSM pipeline: native confocal -> JRC2018U -> JRC2018F
# -> BANC. Adds three Stages on top of the existing chain:
#
#   Stage 0   LSM (Carl Zeiss multi-channel) -> NRRD (NC82 + GFP)
#   Stage A'  native NC82 + GFP -> JRC2018U  (Elastix multi-resolution
#                                             affine + B-spline)
#   Stage B'' JRC2018U -> JRC2018F           (RenderTransformed JAR
#                                             with JRC2018U_JRC2018F.h5;
#                                             dfield is already in the
#                                             output->source direction
#                                             we need, so no `-i`)
#
# Stage C and Stage D are reused from the Kondo pipeline.
# ----------------------------------------------------------------------


#' @title Extract NC82 + GFP channels from a Deng-2019 .lsm file
#' @description Wraps the bundled Python helper
#' \code{inst/python/lsm_to_nrrd.py} (uses \pkg{tifffile} + \pkg{SimpleITK})
#' to extract channel 0 (GFP) and channel 1 (NC82) into two
#' \code{.nrrd} files with proper voxdims metadata. Channel order is
#' assumed; for any LSM cohort that swaps the order, pass the channel
#' indices via \code{...} (forwarded as \code{--gfp-channel} /
#' \code{--nc82-channel}).
#' @param input path to the source .lsm file.
#' @param output_dir directory where \code{<stem>_nc82.nrrd} and
#'   \code{<stem>_gfp.nrrd} are written; created if missing.
#' @param prefix optional output stem (default = input file basename).
#' @param python path to a Python interpreter with \pkg{tifffile} +
#'   \pkg{SimpleITK} installed; defaults to whichever \code{reticulate}
#'   would use.
#' @return a length-2 named character of \code{c(nc82 = ..., gfp = ...)}
#'   absolute paths to the written NRRDs.
#' @export
#' @keywords internal
lsm_to_nrrd <- function(input, output_dir, prefix = NULL, python = NULL,
                        verbose = TRUE) {
  input      <- path.expand(input)
  output_dir <- path.expand(output_dir)
  stopifnot(file.exists(input), nzchar(output_dir))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  helper <- system.file("python", "lsm_to_nrrd.py", package = "neuronbridger")
  if (!nzchar(helper) || !file.exists(helper))
    stop("inst/python/lsm_to_nrrd.py not found in installed package")

  if (is.null(python)) {
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("reticulate is required to find a Python interpreter.")
    python <- reticulate::py_config()$python
  }

  stem <- prefix %||% sub("\\.[^.]+$", "", basename(input))
  args <- c(shQuote(helper), shQuote(input), shQuote(output_dir),
            "--prefix", shQuote(stem))

  rval <- system2(python, args = args,
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("lsm_to_nrrd.py exited ", rval, " on ", basename(input))

  c(nc82 = file.path(output_dir, paste0(stem, "_nc82.nrrd")),
    gfp  = file.path(output_dir, paste0(stem, "_gfp.nrrd")))
}


#' @title Register a native LSM volume to JRC2018U via Elastix
#' @description Runs Elastix multi-resolution affine + B-spline to
#' align a native-confocal NC82 volume to the JRC2018U_HR template,
#' then applies the same transform to the matching GFP channel via
#' \code{transformix}. Parameter files default to the package's
#' \code{inst/extdata/elastix_lm_to_jrc2018u/} bundle.
#' @param nc82 path to native NC82 NRRD (the registration moving image).
#' @param gfp path to native GFP NRRD (same voxdims as \code{nc82}).
#' @param output_dir directory for Elastix outputs (transform params,
#'   warped NC82, warped GFP).
#' @param template path to JRC2018U template NRRD (the registration
#'   fixed image). Defaults to \code{~/templates/JRC2018_UNISEX_20x_HR.nrrd}
#'   --- override via \code{options(neuronbridger.jrc2018u_template = ...)}.
#' @param param_files length-2 character of Elastix parameter file
#'   paths (\code{p_affine.txt}, \code{p_bspline.txt}). Defaults to
#'   the package bundle.
#' @param elastix path to the \code{elastix} binary.
#' @param transformix path to the \code{transformix} binary.
#' @param threads number of threads (Elastix default uses all).
#' @return a list with \code{nc82_jrc2018u} (warped reference NRRD path),
#'   \code{gfp_jrc2018u} (warped GFP NRRD path), \code{params}
#'   (Elastix output directory), and \code{final_metric} (last metric
#'   value reported by Elastix).
#' @export
#' @keywords internal
lm_to_jrc2018u_elastix <- function(nc82,
                                   gfp,
                                   output_dir,
                                   template     = NULL,
                                   param_files  = NULL,
                                   elastix      = "elastix",
                                   transformix  = "transformix",
                                   threads      = 8L,
                                   verbose      = TRUE) {
  nc82       <- path.expand(nc82)
  gfp        <- path.expand(gfp)
  output_dir <- path.expand(output_dir)
  stopifnot(file.exists(nc82), file.exists(gfp), nzchar(output_dir))
  if (Sys.which(elastix) == "")
    stop("`", elastix, "` not on PATH. Install Elastix 5.x.")
  if (Sys.which(transformix) == "")
    stop("`", transformix, "` not on PATH.")

  if (is.null(template))
    template <- getOption("neuronbridger.jrc2018u_template",
                          default = path.expand("~/templates/JRC2018_UNISEX_20x_HR.nrrd"))
  if (!file.exists(template))
    stop("JRC2018U template not found at: ", template,
         ".\nDownload via: ",
         "curl -fL -o ~/templates/JRC2018_UNISEX_20x_HR.nrrd ",
         "https://janelia-flylight-templates.s3.amazonaws.com/JRC2018_Unisex_20x_HR/JRC2018_UNISEX_20x_HR.nrrd")

  if (is.null(param_files)) {
    pdir <- system.file("extdata", "elastix_lm_to_jrc2018u",
                        package = "neuronbridger")
    if (!nzchar(pdir) || !dir.exists(pdir))
      stop("inst/extdata/elastix_lm_to_jrc2018u not in installed package")
    param_files <- c(file.path(pdir, "p_rigid.txt"),
                     file.path(pdir, "p_affine.txt"),
                     file.path(pdir, "p_bspline.txt"))
  }
  for (p in param_files)
    if (!file.exists(p)) stop("missing Elastix parameter file: ", p)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  t0 <- Sys.time()
  args <- c("-f", shQuote(template),
            "-m", shQuote(nc82),
            "-out", shQuote(output_dir),
            "-threads", as.integer(threads))
  for (p in param_files) args <- c(args, "-p", shQuote(p))

  rval <- system2(elastix, args = args,
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("elastix exited ", rval, " --- see ", output_dir, "/elastix.log")

  # Apply the final transform (TransformParameters.<last>.txt) to GFP.
  tps <- list.files(output_dir,
                    pattern = "^TransformParameters\\.[0-9]+\\.txt$",
                    full.names = TRUE)
  if (!length(tps)) stop("no TransformParameters.<n>.txt in ", output_dir)
  final_tp <- tps[order(as.integer(sub("^.*TransformParameters\\.([0-9]+)\\.txt$",
                                       "\\1", tps)))][[length(tps)]]

  gfp_out_dir <- file.path(output_dir, "gfp_xform")
  dir.create(gfp_out_dir, showWarnings = FALSE, recursive = TRUE)
  rval <- system2(transformix,
                  args = c("-in",  shQuote(gfp),
                           "-tp",  shQuote(final_tp),
                           "-out", shQuote(gfp_out_dir),
                           "-threads", as.integer(threads)),
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("transformix (apply-to-GFP) exited ", rval, " --- see ",
         gfp_out_dir, "/transformix.log")

  # Locate result.* (varies by extension among nrrd / mhd)
  nc82_warped <- list.files(output_dir,
                            pattern = "^result\\.[0-9]+\\.(nrrd|mhd|nii(\\.gz)?)$",
                            full.names = TRUE)
  gfp_warped  <- list.files(gfp_out_dir,
                            pattern = "^result\\.(nrrd|mhd|nii(\\.gz)?)$",
                            full.names = TRUE)

  # Pull final metric value out of elastix.log if present
  log_path <- file.path(output_dir, "elastix.log")
  fm <- NA_real_
  if (file.exists(log_path)) {
    tail <- tail(readLines(log_path, warn = FALSE), 200)
    m <- regmatches(tail, regexpr("Final metric value\\s*=\\s*[-+0-9.eE]+", tail))
    if (length(m)) fm <- as.numeric(sub(".*=\\s*", "", m[[length(m)]]))
  }

  if (verbose)
    message(sprintf("  [stage A'] elastix done (%.1fs, final metric = %s)",
                    as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    if (is.na(fm)) "?" else format(fm, digits = 4)))

  list(nc82_jrc2018u = if (length(nc82_warped)) nc82_warped[[length(nc82_warped)]] else NA_character_,
       gfp_jrc2018u  = if (length(gfp_warped))  gfp_warped[[1L]] else NA_character_,
       params        = output_dir,
       final_metric  = fm)
}


#' @title JRC2018U -> JRC2018F image-mode warp via H5 displacement field
#' @description Same shape as \code{\link{fcwb_to_jrc2018f_h5}}, but
#' uses the \code{JRC2018U_JRC2018F.h5} bridging registration. The
#' dfield in that file maps JRC2018F -> JRC2018U; for image
#' resampling we want OUTPUT (JRC2018F) -> SOURCE (JRC2018U), which
#' is exactly the dfield direction --- so unlike \code{fcwb_to_jrc2018f_h5}
#' we DO NOT pass \code{-i} to the JAR.
#' @inheritParams fcwb_to_jrc2018f_h5
#' @export
#' @keywords internal
jrc2018u_to_jrc2018f_h5 <- function(input,
                                    output,
                                    h5_path = NULL,
                                    threads = 8L,
                                    verbose = TRUE) {
  input  <- path.expand(input)
  output <- path.expand(output)
  stopifnot(file.exists(input), nzchar(output))
  jar  <- .neuronbridger_xformimage_jar()
  java <- .neuronbridger_java_bin()
  if (is.null(h5_path)) {
    candidates <- c(
      "~/Library/Application Support/R/nat.jrcbrains/JRC2018U_JRC2018F/JRC2018U_JRC2018F.h5",
      "~/.local/share/R/nat.jrcbrains/JRC2018U_JRC2018F/JRC2018U_JRC2018F.h5"
    )
    h5_path <- Find(file.exists, sapply(candidates, path.expand))
    if (is.null(h5_path))
      stop("JRC2018U_JRC2018F.h5 not found. Run ",
           "nat.jrcbrains::download_saalfeldlab_registrations().")
  }

  interval <- "0,0,0:1651,767,478"        # JRC2018F output grid
  res      <- "0.38,0.38,0.38"
  xfm_arg  <- paste0(h5_path, ":0/dfield:0/invdfield")

  t0 <- Sys.time()
  rval <- system2(java,
                  args = c("-Xmx16g",
                           "-cp", shQuote(jar),
                           "process.RenderTransformed",
                           shQuote(input),
                           shQuote(output),
                           shQuote(interval),
                           "-r", res,
                           "-q", as.integer(threads),
                           shQuote(xfm_arg)),     # NB: no -i here
                  stdout = if (verbose) "" else NULL,
                  stderr = if (verbose) "" else NULL)
  if (rval != 0L)
    stop("RenderTransformed exited ", rval, " on ", basename(input))
  if (!file.exists(output))
    stop("RenderTransformed did not write ", output)

  if (verbose)
    message(sprintf("  [stage B''] %s (%.1fs)",
                    basename(output),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(output)
}


#' @title Bridge a Deng et al. 2019 LSM volume into BANC voxel coordinates
#' @description Top-level orchestrator for the Deng pipeline. Chains
#' Stage 0 (LSM extract) + Stage A' (Elastix native -> JRC2018U) +
#' Stage B'' (H5 JRC2018U -> JRC2018F) + the existing Stage C
#' (transformix JRC2018F -> BANC) + Stage D (nrrd_to_precomputed),
#' producing a Neuroglancer precomputed image layer named
#' \code{<gene>_<sample>.ng/} on the BANC grid.
#'
#' @section Software prerequisites:
#' \itemize{
#'   \item Elastix 5.x \code{elastix} + \code{transformix} on PATH
#'   \item Java 17+ + the saalfeldlab \code{RenderTransformed} JAR
#'   \item Python (via \pkg{reticulate}) with \pkg{tifffile} +
#'         \pkg{SimpleITK}
#'   \item \pkg{nat.jrcbrains} H5 registrations downloaded
#'         (\code{nat.jrcbrains::download_saalfeldlab_registrations()})
#'   \item JRC2018U template NRRD at
#'         \code{~/templates/JRC2018_UNISEX_20x_HR.nrrd}
#'         (download from \code{janelia-flylight-templates} S3).
#' }
#'
#' @param input source .lsm path.
#' @param gene short gene name; forms layer name \code{<gene>_<sample>}.
#' @param sample sample identifier (e.g. \code{BRAIN-F}, \code{VNC-M},
#'   or any user-chosen short tag).
#' @param channel the upstream channel string recorded in the registry
#'   (default \code{"GFP"}).
#' @param dataset master folder under \code{light_level/} (default
#'   \code{"deng_et_al_2019"}).
#' @param output_dir parent directory for intermediates + final
#'   precomputed layer.
#' @param source_path optional original storage path for provenance.
#' @param chunk_size length-3 integer chunk size for the precomputed
#'   layer.
#' @param keep_intermediates logical; keep the per-stage intermediate
#'   NRRDs (default \code{FALSE}).
#' @param threads passed to Elastix + RenderTransformed.
#' @return list with \code{precomputed_dir}, \code{registry_entry}
#'   (one-row tibble), \code{timings} (named numeric, s).
#' @export
#' @keywords internal
lsm_to_banc_layer <- function(input,
                              gene,
                              sample,
                              channel             = "GFP",
                              dataset             = "deng_et_al_2019",
                              output_dir,
                              source_path         = NULL,
                              chunk_size          = c(64L, 64L, 64L),
                              keep_intermediates  = FALSE,
                              threads             = 8L,
                              verbose             = TRUE) {
  stopifnot(file.exists(input), nzchar(gene), nzchar(sample), nzchar(dataset))
  if (missing(output_dir) || is.null(output_dir))
    stop("`output_dir` is required.")
  output_dir <- path.expand(output_dir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  name      <- paste(gene, sample, sep = "_")
  raw_dir   <- file.path(output_dir, "raw")
  reg_dir   <- file.path(output_dir, "elastix", name)
  jrcf_nrrd <- file.path(output_dir, paste0(name, "_in_JRC2018F.nrrd"))
  banc_nrrd <- file.path(output_dir, paste0(name, "_aligned240721_to_BANC.nrrd"))
  pc_dir    <- file.path(output_dir, paste0(name, ".ng"))

  t_all <- Sys.time()
  t0_extract <- system.time({
    chans <- lsm_to_nrrd(input, raw_dir, prefix = name, verbose = verbose)
  })[["elapsed"]]

  t1_register <- system.time({
    r <- lm_to_jrc2018u_elastix(
      nc82       = chans[["nc82"]],
      gfp        = chans[["gfp"]],
      output_dir = reg_dir,
      threads    = threads,
      verbose    = verbose)
  })[["elapsed"]]

  t2_h5 <- system.time(
    jrc2018u_to_jrc2018f_h5(r$gfp_jrc2018u, jrcf_nrrd,
                            threads = threads,
                            verbose = verbose)
  )[["elapsed"]]

  t3_banc <- system.time(
    jrc2018f_to_banc_elastix(jrcf_nrrd, banc_nrrd, verbose = verbose)
  )[["elapsed"]]

  t4_pc <- system.time(
    nrrd_to_precomputed(
      input      = banc_nrrd,
      output     = pc_dir,
      resolution = c(400, 400, 400),
      data_type  = "uint8",
      encoding   = "raw",
      chunk_size = chunk_size,
      overwrite  = TRUE
    )
  )[["elapsed"]]

  if (!keep_intermediates) {
    for (f in c(jrcf_nrrd, banc_nrrd))
      if (file.exists(f)) file.remove(f)
    unlink(reg_dir, recursive = TRUE)
    unlink(raw_dir, recursive = TRUE)
  }

  registry_entry <- tibble::tibble(
    dataset         = dataset,
    gene            = gene,
    sample          = sample,
    channel         = channel,
    name            = name,
    gs_url          = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/%s/%s.ng/",
                              dataset, name),
    alignment_space = "BANC",
    voxdims_nm      = list(c(400L, 400L, 400L)),
    source_filename = basename(input),
    source_path     = source_path %||% input,
    uploaded        = format(Sys.Date())
  )

  list(precomputed_dir = pc_dir,
       registry_entry  = registry_entry,
       elastix_metric  = r$final_metric,
       timings         = c(stage0_lsm_extract       = t0_extract,
                           stageA_native_to_jrcu    = t1_register,
                           stageB_jrcu_to_jrcf      = t2_h5,
                           stageC_jrcf_to_banc      = t3_banc,
                           stageD_precomputed       = t4_pc,
                           total                    = as.numeric(
                             difftime(Sys.time(), t_all, units = "secs"))))
}
