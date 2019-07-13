

context("Generate exchange partners")
library("anticlust")

test_that("Exchange partners are generated correctly", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- NULL
  preclusters <- NULL
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(all(partners == c(2, 4, 6)), TRUE)

  ## Add a categorical restriction
  categories <- c(1, 1, 1, 1, 2, 2)
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(all(partners == c(2, 4)), TRUE)

  ## Add a preclustering restriction
  categories <- NULL
  preclusters <- c(1, 1, 2, 2, 3, 3)
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(partners == 2, TRUE)

  ## Use both restrictions
  categories <- c(1, 1, 1, 1, 2, 2)
  preclusters <- c(1, 2, 2, 1, 3, 3)
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(partners == 4, TRUE)
})

test_that("Exchange partners are generated correctly when NAs are present", {
  clusters <- c(1, 2, 1, 2, 1, 2, NA, NA)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- NULL
  preclusters <- NULL
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(all(partners == c(2, 4, 6, 7, 8)), TRUE)

  ## NAs and preclustering restrictions
  ## (a) NA elements not in the same precluster
  preclusters <- c(1, 1, 2, 2, 3, 3, 4, 4)
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(partners == 2, TRUE)
  ## (b) One NA element in the same precluster
  preclusters <- c(1, 1, 2, 2, 3, 3, 1, 4)
  partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  expect_equal(all(partners == c(2 ,7)), TRUE)

  ## NAs and categorical restrictions - TODO
  ##preclusters <- NULL
  ##categories <- c(1, 1, 1, 1, 2, 2, 2, 1)
  ##partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  ##expect_equal(partners == c(2, 4, 8), TRUE)
  ## (b) One NA element in the same precluster
  ##preclusters <- c(1, 1, 2, 2, 3, 3, 1, 4)
  ##partners <- get_exchange_partners(clusters, i, group_i, N, categories, preclusters)
  ##expect_equal(all(partners == c(2 ,7)), TRUE)
})
