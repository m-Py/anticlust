library("anticlust")

context("New C implementation")

test_that("New k-means C implementation yields same results as other implementations", {
  
  set.seed(123)
  M <- 5
  N <- 180
  K <- 4
  features <- matrix(rnorm(N*M), ncol = M)
  init <- sample(rep_len(1:K, nrow(features)))
  

  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = Inf)
  ac_exchangeC2 <- anticlustering(features, K = init, objective = "variance")

  expect_true(all(ac_exchangeC == ac_exchangeC2))
  
  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = 20)
  ac_exchangeR <- fast_anticlustering(features, K = init, backend = "R", k_neighbours = 20)

  expect_true(all(ac_exchangeC == ac_exchangeR))
  
  # Also test with categorical restrictions
  features <- schaper2019[, 3:6]
  categories <- schaper2019$room
  init <- initialize_clusters(nrow(features), K = 3, categories)
  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = Inf, categories = categories)
  ac_exchangeC2 <- anticlustering(features, K = init, objective = "variance", categories = categories)
  expect_true(all(ac_exchangeC == ac_exchangeC2))
  
  # Use reduced exchange partners
  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = 10, categories = categories)
  ac_exchangeR <- fast_anticlustering(features, K = init, backend = "R", k_neighbours = 10, categories = categories)
  expect_true(all(ac_exchangeC == ac_exchangeR))
  
  # What if `k_neighbours` is greater than the number of elements in the group with fewest members?
  categories <- merge_into_one_variable(cbind(schaper2019$syllables, schaper2019$room))
  init <- initialize_clusters(nrow(features), K = 3, categories)
  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = 10, categories = categories)
  ac_exchangeR <- fast_anticlustering(features, K = init, backend = "R", k_neighbours = 10, categories = categories)
  expect_true(all(ac_exchangeC == ac_exchangeR))
  
})
