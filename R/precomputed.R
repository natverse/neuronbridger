#' @title Convert an .nrrd image volume to Neuroglancer "precomputed" format
#'
#' @description Reads a 3D image volume from an \code{.nrrd} file (or accepts
#' one already loaded into R) and writes it as a Neuroglancer
#' \href{https://github.com/google/neuroglancer/tree/master/src/neuroglancer/datasource/precomputed}{precomputed}
#' image layer — an \code{info} JSON plus chunked, optionally-gzipped raw
#' bricks under \code{<resolution>/<x_min-x_max>_<y_min-y_max>_<z_min-z_max>}.
#' The output can be opened directly by \href{https://spelunker.cave-explorer.org/}{spelunker}
#' or any other Neuroglancer viewer when served over HTTP/S3/GCS, and is the
#' format the BANC connectome and most CAVE-served layers use.
#'
#' This is a thin R wrapper around
#' \href{https://github.com/seung-lab/cloud-volume}{cloud-volume}: each
#' \code{nrrd_to_precomputed()} call boots a single Python process via
#' \pkg{reticulate}, builds an info dict, and writes the volume as a single
#' bounding-box write. For very large volumes you'll usually want to chain
#' \code{cloudvolume.transfer.TransferTask} or generate downsampled mip
#' levels with \code{cloudvolume.tasks.create_downsampling_tasks}; both are
#' available from the same Python session via \code{reticulate::import("cloudvolume")}.
#'
#' Users with an existing Python pipeline may prefer the
#' \href{https://github.com/jasper-tms/npimage/blob/main/npimage/imageio.py}{\code{npimage.save}}
#' one-liner, which wraps the same \code{cloud-volume} call:
#' \code{npimage.save(arr, 'foo.ng', pixel_size = res_nm)} for
#' single-channel \code{.nrrd}-style inputs. Multi-channel light-sheet
#' \code{.lsm} stacks need a per-channel loop or RGB packing in either
#' approach.
#'
#' @param input either a path to an \code{.nrrd} file or a 3D numeric/integer
#' array already loaded into R. If a path, voxel size is read from the NRRD
#' header (\code{space directions}); for an in-memory array, supply
#' \code{resolution} explicitly.
#' @param output where to write the precomputed layer. A local directory
#' path is accepted (\code{file://} prefix is added if missing); pass a
#' \code{gs://}, \code{s3://} or \code{https://} URL to write directly to a
#' bucket (requires the relevant cloud-volume credentials).
#' @param resolution numeric vector of length 3, voxel size in
#' \strong{nanometres}. If \code{NULL} (default) and \code{input} is a
#' path, voxel size is taken from the NRRD header (assumed to be in microns
#' and multiplied by 1000). If the array has been down-sampled before
#' calling, multiply the source resolution by the down-sampling factor.
#' @param data_type one of \code{"uint8"} (default), \code{"uint16"},
#' \code{"uint32"}, \code{"float32"}. Choose to match the dynamic range of
#' your data: 16-bit microscopy data with most signal in the bottom 12 bits
#' typically fits cleanly into \code{"uint8"} after a right-shift of 4 bits.
#' @param encoding one of \code{"raw"} (default; gzip'd at the chunk
#' boundary), \code{"compresso"}, \code{"jpeg"} (lossy; image layers only),
#' or \code{"compressed_segmentation"} (segmentation layers only).
#' @param chunk_size length-3 integer vector of voxel chunk dimensions.
#' Default \code{c(64, 64, 32)}.
#' @param voxel_offset length-3 integer vector applied as the global voxel
#' origin in the Neuroglancer view. Default \code{c(0, 0, 0)}.
#' @param layer_type one of \code{"image"} (default) or \code{"segmentation"}.
#' @param compress logical; gzip compression for raw chunks (default
#' \code{TRUE}). Ignored for \code{encoding != "raw"}.
#' @param overwrite logical; if \code{FALSE} (default) and \code{output} is
#' a local directory that already exists, abort.
#' @param python optional path or virtualenv identifier passed to
#' \code{reticulate::use_python()} before importing \code{cloudvolume}.
#'
#' @return Invisibly, the resolved \code{output} URL (with \code{file://}
#' prefix added if needed).
#'
#' @examples
#' \dontrun{
#' # Convert a Kondo-2020 receptor expression NRRD (IS2-space, 16-bit) into a
#' # Neuroglancer precomputed layer down-sampled 4x in xy and 2x in z, mapped
#' # to 8-bit:
#' library(nat); library(reticulate); py_require("cloud-volume")
#'
#' v <- read.nrrd("IS2_CapaR_no1_01_warp_m0g40c4e1e-1x16r3.nrrd")
#' hdr <- attr(v, "header")
#' src_res_um <- diag(hdr[["space directions"]])
#'
#' v8 <- v[seq(1, dim(v)[1], by = 4),
#'         seq(1, dim(v)[2], by = 4),
#'         seq(1, dim(v)[3], by = 2)]
#' v8 <- as.integer(pmin(pmax(v8, 0), 4095) / 16)  # 12-bit -> 8-bit
#' dim(v8) <- c(192, 192, 87)
#'
#' nrrd_to_precomputed(
#'   v8,
#'   output     = "file:///tmp/capar_pc",
#'   resolution = round(src_res_um * c(4, 4, 2) * 1000),
#'   data_type  = "uint8"
#' )
#' # Then upload:
#' #   gsutil cp -r /tmp/capar_pc gs://your-bucket/lm_layers/CapaR_IS2/
#' }
#' @seealso \code{bancr::banc_lm_scene} for assembling a BANC-overlay
#' Neuroglancer state pointing at a precomputed layer; \code{\link{nrrd_to_mip}}.
#' @export
nrrd_to_precomputed <- function(input,
                                output,
                                resolution   = NULL,
                                data_type    = c("uint8", "uint16",
                                                 "uint32", "float32"),
                                encoding     = c("raw", "compresso",
                                                 "jpeg",
                                                 "compressed_segmentation"),
                                chunk_size   = c(64L, 64L, 32L),
                                voxel_offset = c(0L, 0L, 0L),
                                layer_type   = c("image", "segmentation"),
                                compress     = TRUE,
                                overwrite    = FALSE,
                                python       = NULL) {
  data_type  <- match.arg(data_type)
  encoding   <- match.arg(encoding)
  layer_type <- match.arg(layer_type)
  if (length(chunk_size)   != 3L) stop("`chunk_size` must be length 3.")
  if (length(voxel_offset) != 3L) stop("`voxel_offset` must be length 3.")

  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Install the 'reticulate' package: install.packages('reticulate').")
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)

  # Resolve input
  if (is.character(input) && length(input) == 1L && file.exists(input)) {
    if (!grepl("\\.nrrd$", input, ignore.case = TRUE))
      stop("Expected an .nrrd file: ", input)
    if (!requireNamespace("nat", quietly = TRUE))
      stop("Install the 'nat' package to read .nrrd input.")
    vol <- nat::read.nrrd(input)
    if (is.null(resolution)) {
      sd <- attr(vol, "header")[["space directions"]]
      resolution <- as.integer(round(diag(sd) * 1000))  # um -> nm
    }
  } else if (is.array(input) || is.numeric(input)) {
    vol <- input
    if (is.null(resolution))
      stop("`resolution` (nm) must be supplied when `input` is an in-memory array.")
  } else {
    stop("`input` must be a path to an .nrrd file or a 3D numeric array.")
  }
  if (length(dim(vol)) != 3L)
    stop("Expected a 3D volume; got dim = ", paste(dim(vol), collapse = " x "))
  if (length(resolution) != 3L)
    stop("`resolution` must be length 3 (x, y, z) in nanometres.")

  # Resolve output URL: bare local path -> file://
  if (!grepl("^[a-z0-9+.-]+://", output, ignore.case = TRUE))
    output <- paste0("file://", normalizePath(output, mustWork = FALSE))
  is_local <- startsWith(output, "file://")
  if (is_local) {
    local_path <- substring(output, 8L)
    if (dir.exists(local_path)) {
      if (!overwrite)
        stop("`output` directory exists: ", local_path,
             " (pass overwrite = TRUE to replace it)")
      unlink(local_path, recursive = TRUE)
    }
    dir.create(local_path, showWarnings = FALSE, recursive = TRUE)
  }

  # Cast to the requested numeric type
  storage_mode <- switch(data_type,
                         uint8 = "integer", uint16 = "integer",
                         uint32 = "integer", float32 = "double")
  storage.mode(vol) <- storage_mode

  np <- reticulate::import("numpy",       convert = FALSE)
  cv <- reticulate::import("cloudvolume", convert = FALSE)

  arr <- np$asarray(reticulate::r_to_py(vol), dtype = data_type)

  info <- cv$CloudVolume$create_new_info(
    num_channels = 1L,
    layer_type   = layer_type,
    data_type    = data_type,
    encoding     = encoding,
    resolution   = as.list(as.integer(round(resolution))),
    voxel_offset = as.list(as.integer(voxel_offset)),
    chunk_size   = as.list(as.integer(chunk_size)),
    volume_size  = as.list(as.integer(dim(vol)))
  )

  vol_obj <- cv$CloudVolume(output, info = info, mip = 0L,
                            compress = if (compress && encoding == "raw")
                                         "gzip" else FALSE)
  vol_obj$commit_info()
  bb <- cv$Bbox$from_list(c(0L, 0L, 0L, dim(vol)[1], dim(vol)[2], dim(vol)[3]))
  vol_obj[bb] <- arr

  message("Wrote precomputed layer to: ", output)
  invisible(output)
}

# BANC-specific scene assembly (`banc_lm_scene()`) lives in the `bancr`
# package, which exists for that purpose. See
# vignettes/lm_layer_neuroglancer.Rmd for the worked example.
