skip_if_not_installed("reticulate")
skip_if_not(reticulate::py_module_available("cloudvolume"),
            "cloudvolume Python package not available")

test_that("nrrd_to_precomputed round-trips a uint8 volume through cloud-volume", {
  td <- tempfile(); on.exit(unlink(td, recursive = TRUE), add = TRUE)

  set.seed(11)
  v <- array(as.integer(runif(64 * 32 * 16, 0, 200)),
             dim = c(64L, 32L, 16L))

  out <- nrrd_to_precomputed(v,
                             output     = td,
                             resolution = c(1000, 1000, 2000),
                             data_type  = "uint8",
                             chunk_size = c(32L, 32L, 16L))
  expect_true(startsWith(out, "file://"))
  expect_true(file.exists(file.path(td, "info")))

  info <- jsonlite::fromJSON(file.path(td, "info"))
  expect_equal(info$type,      "image")
  expect_equal(info$data_type, "uint8")
  expect_equal(info$scales$resolution[[1]], c(1000, 1000, 2000))
  expect_equal(info$scales$size[[1]],       c(64, 32, 16))

  np  <- reticulate::import("numpy",       convert = TRUE)
  cv  <- reticulate::import("cloudvolume", convert = FALSE)
  vol <- cv$CloudVolume(out, mip = 0L, fill_missing = TRUE)
  arr <- np$squeeze(np$asarray(vol[0:64, 0:32, 0:16]), axis = 3L)
  expect_equal(as.integer(arr), as.integer(v))
})

# Note: BANC-specific scene assembly (banc_lm_scene) lives in the bancr
# package and is tested there.
