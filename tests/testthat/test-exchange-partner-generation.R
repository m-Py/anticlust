

context("Generate exchange partners")
library("anticlust")

test_that("Exchange partners are generated correctly - no restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- NULL
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4, 6)), TRUE)
})

test_that("Exchange partners are generated correctly - categorical restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- c(1, 1, 1, 1, 2, 2)
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4)), TRUE)
})

test_that("Exchange partners are generated correctly - preclustering restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  # Add a preclustering restriction (this is now just the same as
  # categorical restrictions and does not really make sense here)
  preclusters <- c(1, 1, 2, 2, 3, 3)
  partners <- get_exchange_partners(clusters, i, preclusters)
  expect_equal(partners == 2, TRUE)
})

test_that("Exchange partners are generated correctly - preclustering and categorical restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  ## Use both restrictions
  categories <- c(1, 1, 1, 1, 2, 2)
  preclusters <- c(1, 2, 2, 1, 3, 3)
  constraints <- merge_into_one_variable(cbind(categories, preclusters))
  partners <- get_exchange_partners(clusters, i, constraints)
  expect_equal(partners == 4, TRUE)
})

test_that("Exchange partners are generated correctly when NAs are present - no restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2, NA, NA)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- NULL
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4, 6, 7, 8)), TRUE)
})

test_that("Exchange partners are generated correctly when NAs are present - preclustering restruction", {
  clusters <- c(1, 2, 1, 2, 1, 2, NA, NA)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- NULL
  ## NAs and preclustering restrictions
  ## (a) NA elements not in the same precluster
  preclusters <- c(1, 1, 2, 2, 3, 3, 4, 4)
  partners <- get_exchange_partners(clusters, i, preclusters)
  expect_equal(partners == 2, TRUE)
  ## (b) One NA element in the same precluster
  preclusters <- c(1, 1, 2, 2, 3, 3, 1, 4)
  partners <- get_exchange_partners(clusters, i, preclusters)
  expect_equal(all(partners == c(2 ,7)), TRUE)
})

test_that("Exchange partners are generated correctly when NAs are present - categorical restruction", {
  clusters <- c(1, 2, 1, 2, 1, 2, NA, NA)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  preclusters <- NULL
  ## (a) No NA element in the same category
  categories <- c(1, 1, 1, 1, 2, 2, 2, 2)
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4)), TRUE)
  ## (b) One NA elements is in the same category
  categories <- c(1, 1, 1, 1, 2, 2, 2, 1)
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4, 8)), TRUE)
})
