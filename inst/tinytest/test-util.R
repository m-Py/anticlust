library("anticlust")

context("Utility function tests")

test_that("test of input works as expected", {
  features <- matrix(rnorm(12), nrow = 4)
  ## False number of anticlusters:
  expect_error(legal_number_of_clusters(features, 1))
  ## Only one distinct anticluster:
  expect_error(legal_number_of_clusters(features, rep(1, nrow(features))))
  ## Anticlusters do not occur equally often:
  expect_error(legal_number_of_clusters(features, c(1, 1, 1, 2)))
  correct <- legal_number_of_clusters(matrix(rnorm(12), nrow = 4), c(1, 2, 1, 2))
  expect_equal(correct, NULL)
})
