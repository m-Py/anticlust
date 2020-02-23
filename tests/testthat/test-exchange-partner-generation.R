

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


test_that("Exchange partners are generatede correctly for fast k-means method", {
  # function that tests that all exchange partners are from the same category
  test_idx <- function(idx, categories) {
    if (argument_exists(categories)) {
      for (i in 1:length(idx)) {
        partners <- categories[idx[[i]]]
        expect_true(all(partners == partners[1]))
      }
    }
  }
  
  ### Using categorical restrictions
  # Case 1: restricted number of exchange partners (i.e., nearest neighbour search)
  categories <- to_numeric(schaper2019$room)
  partners <- all_exchange_partners(
    features = schaper2019[, 3:6],
    k_neighbours = 10,
    categories = categories
  )
  test_idx(partners, categories)
  ## ensure correct order of the output 
  expect_true(all(sapply(partners, function(x) x[1]) == 1:96))
  
  # Case 2: no restriction on number of exchange partners
  partners <- all_exchange_partners(
    features = schaper2019[, 3:6],
    k_neighbours = Inf,
    categories = categories
  )
  test_idx(partners, categories)
  
  ### Not using categorical restrictions
  # Case 1: restricted number of exchange partners (i.e., nearest neighbour search)
  categories <- to_numeric(schaper2019$room)
  k_neighbours <- 10
  partners <- all_exchange_partners(
    features = schaper2019[, 3:6],
    k_neighbours = k_neighbours,
    categories = NULL
  )
  expect_true(all(sapply(partners, length) == k_neighbours + 1))
  ## ensure correct order of the output 
  expect_true(all(sapply(partners, function(x) x[1]) == 1:96))
  
  # Case 2: no restriction on number of exchange partners
  partners <- all_exchange_partners(
    features = schaper2019[, 3:6],
    k_neighbours = Inf,
    categories = NULL
  )
  expect_true(all(sapply(partners, length) == nrow(schaper2019)))
})
