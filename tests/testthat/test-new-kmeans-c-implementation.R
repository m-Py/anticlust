library("anticlust")

context("New C implementation")

test_that("New k-means C implementation yields same results as other implementations", {

  M <- 5
  N <- 180
  K <- 4
  features <- matrix(rnorm(N*M), ncol = M)
  init <- sample(rep_len(1:K, nrow(features)))
  

  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = Inf)
  ac_exchangeC2 <- anticlustering(features, K = init, objective = "variance")

  expect_true(all(ac_exchangeC == ac_exchangeC2))
  
  ac_exchangeC <- fast_anticlustering(features, K = init, backend = "C", k_neighbours = 20, nn_method = "RANN")
  ac_exchangeR <- fast_anticlustering(features, K = init, backend = "R", k_neighbours = 20, nn_method = "RANN")

  expect_true(all(ac_exchangeC == ac_exchangeR))

})
