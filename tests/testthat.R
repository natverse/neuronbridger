library(testthat)
library(neuronbridger)

conn=try(neuronbridge_ids())
if(isTRUE(length(conn))) {
  Sys.setenv(SKIP_NP_SERVER_TESTS="FALSE")
} else {
  Sys.setenv(SKIP_NP_SERVER_TESTS="TRUE")
}

test_check("neuronbridger")
