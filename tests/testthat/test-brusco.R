
context("Brusco et al. algorithm")
library("anticlust")

test_that("input and output work expectedly for Brusco algorithm", {
  
  # Test standard usage: K = number of groups
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = 3,
    method = "brusco"
  )
  expect_true(all(table(anticlusters) == c(32, 32, 32)))
  
  # Different group sizes 
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = c(48, 24, 24),
    method = "brusco"
  )
  expect_true(all(table(anticlusters) == c(48, 24, 24)))

  # no categorical restrictions
  expect_error(
    anticlusters <- anticlustering(
      schaper2019[, 3:6],
      K = 4,
      categories = schaper2019$room,
      method = "brusco"
    ), 
    regexp = "categorical"
  )
  
  # no categorical restrictions
  expect_error(
    anticlusters <- anticlustering(
      schaper2019[, 3:6],
      K = 4,
      preclustering = TRUE,
      method = "brusco"
    ), 
    regexp = "preclustering"
  )
  
})
