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
#' @param verbose logical; emit stage timing.
#' @return path to the FCWB-aligned NRRD.
#' @export
is2_to_fcwb_cmtk <- function(input,
                             output,
                             sample_space = "IS2",
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

  # FCWB grid: 1769 x 1026 x 108 at 0.318967, 0.318427, 1.0 um.
  # CMTK's parser balks on >7 decimal-place voxel sizes, so 6 dp.
  target_grid <- "1769,1026,108:0.318967,0.318427,1.0"

  t0 <- Sys.time()
  rval <- system2(reformatx,
                  args = c("--target-grid", shQuote(target_grid),
                           "--floating",    shQuote(input),
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
                           shQuote(xfm_arg)),
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
#' @param sample sample identifier (e.g. \code{"01"}).
#' @param channel source-stack channel string (default \code{"no2"}).
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
#' @param threads passed to \code{fcwb_to_jrc2018f_h5}.
#' @return list with \code{precomputed_dir} (filesystem path),
#'   \code{registry_entry} (one-row tibble), and \code{timings}
#'   (named numeric, seconds per stage).
#' @export
lm_to_banc_layer <- function(input,
                             gene,
                             sample,
                             channel             = "no2",
                             dataset             = "kondo_et_al_2020",
                             output_dir,
                             source_path         = NULL,
                             chunk_size          = c(64L, 64L, 64L),
                             keep_intermediates  = FALSE,
                             sample_space        = "IS2",
                             threads             = 8L,
                             verbose             = TRUE) {
  stopifnot(file.exists(input),
            nzchar(gene), nzchar(sample), nzchar(dataset))
  if (missing(output_dir) || is.null(output_dir))
    stop("`output_dir` is required.")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  name      <- paste(gene, sample, sep = "_")
  fcwb_nrrd <- file.path(output_dir, paste0(name, "_in_FCWB.nrrd"))
  jrcf_nrrd <- file.path(output_dir, paste0(name, "_in_JRC2018F.nrrd"))
  banc_nrrd <- file.path(output_dir, paste0(name, "_aligned240721_to_BANC.nrrd"))
  pc_dir    <- file.path(output_dir, paste0(name, ".ng"))

  t_all <- Sys.time()
  tA <- system.time(is2_to_fcwb_cmtk(input, fcwb_nrrd,
                                     sample_space = sample_space,
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
