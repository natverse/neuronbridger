test_that("colormip_from_array reproduces the BANC depth_lut algorithm", {
  # Three points at known z-depths in a JRC2018U_HR-shaped volume
  vol <- array(0L, dim = c(20, 20, 174))
  vol[2,  5, 1]   <- 255L  # most anterior
  vol[5,  5, 87]  <- 255L  # middle (z=87 -> int((86/173)*255)=126 -> LUT[127])
  vol[8,  5, 174] <- 255L  # most posterior

  cmip <- colormip_from_array(vol, target_space = "brain")

  # Output orientation: rows = y, cols = x, channels = 3
  expect_equal(dim(cmip), c(20, 20, 3))

  # Pixel-perfect match to fanc.render_neurons.depth_lut
  expect_equal(round(cmip[5, 2, ] * 255), c(127, 0, 255))   # LUT[1]
  expect_equal(round(cmip[5, 5, ] * 255), c(123, 255, 132)) # LUT[127]
  expect_equal(round(cmip[5, 8, ] * 255), c(255, 0, 200))   # LUT[256]

  # Background remains black
  expect_equal(round(cmip[1, 1, ] * 255), c(0, 0, 0))
})

test_that("colormip_from_array prepends a 90-pixel header in VNC space", {
  vol <- array(0L, dim = c(30, 30, 219))
  vol[10, 10, 50] <- 255L
  cmip <- colormip_from_array(vol, target_space = "VNC")
  expect_equal(dim(cmip), c(30 + 90, 30, 3))
  expect_true(all(cmip[1:90, , ] == 0))  # header is black
})

test_that("nrrd_to_mip round-trips through PNG and dispatches on folders", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  v <- array(0L, dim = c(40, 40, 174)); v[20, 20, 87] <- 255L
  nat::write.nrrd(v, file.path(td, "test.nrrd"))

  out <- nrrd_to_mip(td, target_space = "brain")
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))

  img <- png::readPNG(out[[1]])
  expect_equal(dim(img), c(40, 40, 3))
  expect_equal(round(img[20, 20, ] * 255), c(123, 255, 132))
})
