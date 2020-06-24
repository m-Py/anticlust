
context("Test exchange methods")
library("anticlust")

test_that("fast exchange and exchange functions yield the same results - kmeans", {
  # for K = 2
  clusters <- rep(1:2, 5)
  features <- matrix(rnorm(20), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = variance_objective)
  ac_fast <- fast_anticlustering(features, clusters)
  expect_equal(all(ac == ac_fast), TRUE)

  ## for K = 3
  clusters <- rep(1:3, 5)
  features <- matrix(rnorm(30), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = variance_objective)
  ac_fast <- fast_anticlustering(features, clusters)
  expect_equal(all(ac == ac_fast), TRUE)

  ## for K = 4
  clusters <- rep(1:4, 5)
  features <- matrix(rnorm(20), ncol = 1)
  ac <- anticlustering(features, K = clusters, objective = variance_objective)
  ac_fast <- fast_anticlustering(features, clusters)
  expect_equal(all(ac == ac_fast), TRUE)
})


test_that("fast exchange and exchange functions yield the same results - anticluster editing", {
  # for K = 2
  K <- 2
  clusters <- rep(1:K, 5)
  features <- matrix(rnorm(K * 10), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = diversity_objective)
  ac_fast <- fast_exchange_dist(as.matrix(dist(features)), clusters, NULL)
  expect_equal(all(ac == ac_fast), TRUE)

  # for K = 3
  K <- 3
  clusters <- rep(1:K, 5)
  features <- matrix(rnorm(K * 10), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = diversity_objective)
  ac_fast <- fast_exchange_dist(as.matrix(dist(features)), clusters, NULL)
  expect_equal(all(ac == ac_fast), TRUE)

  # for K = 4
  K <- 4
  clusters <- rep(1:K, 5)
  features <- matrix(rnorm(K * 10), ncol = 2)
  ac <- anticlustering(features, K = clusters, objective = diversity_objective)
  ac_fast <- fast_exchange_dist(as.matrix(dist(features)), clusters, NULL)
  expect_equal(all(ac == ac_fast), TRUE)

  ## Test preclustering restrictions
  features <- schaper2019[, 3:6]
  clusters <- categorical_sampling(matching(features, 2), K = 2)
  ac <- anticlustering(features, K = clusters, objective = diversity_objective)
  ac_fast <- anticlustering(features, K = clusters, objective = "distance")
  expect_equal(all(ac == ac_fast), TRUE)

  ## Test categorical restrictions
  features <- schaper2019[, 3:6]
  ac <- anticlustering(features, K = rep_len(1:2, nrow(schaper2019)),
                       objective = diversity_objective,
                       categories = schaper2019$room)
  ac_fast <- anticlustering(features, K = rep_len(1:2, nrow(schaper2019)),
                            objective = "distance",
                            categories = schaper2019$room)
  expect_equal(all(ac == ac_fast), TRUE)

})
