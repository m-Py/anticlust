

context("Test solving function for distance anticlustering")
library("anticlust")

test_that("all levels of heuristicism work and that exact approach has best objective", {
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    p_anticlusters <- conditions[k, "p"]
    n_elements <- p_anticlusters * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    ## Traverse through levels of heuristicism
    obj_values <- rep(NA, 4)
    anti_list <- list()
    for (i in c(TRUE, FALSE)) {
      anticlusters <- exact_anticlustering(as.matrix(dist(features)),
                                           p_anticlusters, preclustering = i)
      anti_list[[i + 1]] <- anticlusters
      # Allow for some numeric imprecision of ILP solver:
      obj_values[i + 1]  <- round(diversity_objective_(anticlusters, features), 10)
    }
    ## Exact solution must have maximum objective
    expect_equal(which.max(obj_values), 1)
  }
})

test_that("Solving ILP works as expected for max dispersion", {
  N <- 30
  K <- 3
  M <- 5
  dat <- matrix(rnorm(N*M), ncol = M)
  
  ILP <- anticlustering(dat, K = K, objective = "dispersion", method = "ilp")
  HEURISTIC <- anticlustering(dat, K = K, objective = "dispersion", method = "local-maximum")
  
  expect_true(dispersion_objective(dat, ILP) >= dispersion_objective(dat, HEURISTIC))
  expect_true(all(table(ILP) == N/K))
})

