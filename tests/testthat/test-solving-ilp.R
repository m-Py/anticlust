

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
