
context("Test exchange methods")
library("anticlust")

test_that("fast exchange and exchange functions yield the same results", {
  # for K = 2
  clusters <- rep(1:2, 5)
  features <- matrix(rnorm(20), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = variance_objective_)
  ac_fast <- fast_anticlustering(features, clusters)
  expect_equal(all(ac == ac_fast), TRUE)

  ## for K = 3
  clusters <- rep(1:3, 5)
  features <- matrix(rnorm(30), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = variance_objective_)
  ac_fast <- fast_anticlustering(features, clusters)
  expect_equal(all(ac == ac_fast), TRUE)

  ## for K = 4
  clusters <- rep(1:4, 5)
  features <- matrix(rnorm(20), ncol = 1)
  ac <- anticlustering(features, K = clusters, objective = variance_objective_)
  ac_fast <- fast_anticlustering(features, clusters)
  expect_equal(all(ac == ac_fast), TRUE)
})
