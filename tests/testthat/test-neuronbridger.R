skip_if_offline()
skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("neuronbridger works", {

  expect_is(nb.a <- neuronbridge_avoid("542634818","754538881"), 'data.frame')

  expect_is(nb.a <- neuronbridge_predict_split("R84D10","R26B04"), 'data.frame')

  expect_is(nb.mip <- neuronbridge_mip("542634818"), 'array')

  expect_is(scan.mip <- scan_mip("542634818",sleep=0.1,type="id"), NULL)

  unlink(unlist(options("neuronbridger")), recursive = TRUE)

})
