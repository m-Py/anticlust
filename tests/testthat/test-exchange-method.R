
context("Test exchange methods")
library("anticlust")

test_that("fast exchange and exchange functions yield the same results - kmeans", {
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


test_that("fast exchange and exchange functions yield the same results - anticluster editing", {
  # for K = 2
  K <- 2
  clusters <- rep(1:K, 5)
  features <- matrix(rnorm(K * 10), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = obj_value_distance)
  ac_fast <- fast_exchange_dist(as.matrix(dist(features)), clusters, NULL)
  expect_equal(all(ac == ac_fast), TRUE)

  # for K = 3
  K <- 3
  clusters <- rep(1:K, 5)
  features <- matrix(rnorm(K * 10), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = obj_value_distance)
  ac_fast <- fast_exchange_dist(as.matrix(dist(features)), clusters, NULL)
  expect_equal(all(ac == ac_fast), TRUE)

  # for K = 4
  K <- 4
  clusters <- rep(1:K, 5)
  features <- matrix(rnorm(K * 10), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = obj_value_distance)
  ac_fast <- fast_exchange_dist(as.matrix(dist(features)), clusters, NULL)
  expect_equal(all(ac == ac_fast), TRUE)

  ## Test preclustering restrictions
  features <- schaper2019[, 3:6]
  ac <- anticlustering(features, K = rep_len(1:2, nrow(schaper2019)),
                       objective = obj_value_distance, preclustering = TRUE)
  ac_fast <- anticlustering(features, K = rep_len(1:2, nrow(schaper2019)),
                            objective = "distance", preclustering = TRUE)
  expect_equal(all(ac == ac_fast), TRUE)

  ## Test categorical restrictions
  features <- schaper2019[, 3:6]
  ac <- anticlustering(features, K = rep_len(1:2, nrow(schaper2019)),
                       objective = obj_value_distance,
                       categories = schaper2019$room)
  ac_fast <- anticlustering(features, K = rep_len(1:2, nrow(schaper2019)),
                            objective = "distance",
                            categories = schaper2019$room)
  expect_equal(all(ac == ac_fast), TRUE)

  ## Test categorical AND preclustering restrictions [not possible at the moment!]
  features <- schaper2019[, 3:6]
  expect_error(
    anticlustering(
      features,
      K = rep_len(1:2, nrow(schaper2019)),
      objective = obj_value_distance,
      preclustering = TRUE,
      categories = schaper2019$room)
  )
  expect_error(
    anticlustering(
      features,
      K = rep_len(1:2, nrow(schaper2019)),
      objective = "distance",
      preclustering = TRUE,
      categories = schaper2019$room
    )
  )
})
