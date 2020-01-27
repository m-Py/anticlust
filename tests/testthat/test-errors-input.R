
context("Error handling anticlustering")
library("anticlust")

test_that("exported function cannot be used with NA", {
  data(iris)
  iris[1, 1] <- NA
  expect_error(
    anticlustering(iris[, -5], K = 3)
  )
  expect_error(
    matching(iris[, -5], p = 3)
  )
  expect_error(
    balanced_clustering(iris[, -5], p = 3)
  )
})

test_that("exported functions cannot be used with non-numeric input", {
  data(iris)
  expect_error(
    anticlustering(iris, K = 3)
  )
  expect_error(
    matching(iris, p = 3)
  )
  expect_error(
    balanced_clustering(iris, p = 3)
  )
})
