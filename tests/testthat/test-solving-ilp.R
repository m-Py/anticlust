

context("Test solving function for distance anticlustering")
library("anticlust")

test_that("all levels of heuristicism work and that exact approach has best objective", {
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    p_anticlusters <- conditions[k, "p"]
    n_elements <- p_anticlusters * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)
    solver <- "Rglpk"
    ## Traverse through levels of heuristicism
    obj_values <- rep(NA, 4)
    anti_list <- list()
    for (i in 0:1) {
      anticlusters <- distance_anticlustering(features, p_anticlusters,
                                              solver, heuristic = i)
      anti_list[[i + 1]] <- anticlusters
      # Allow for some numeric imprecision of ILP solver:
      obj_values[i + 1]  <- round(get_objective(features, anticlusters, "distance"), 10)
    }
    ## Exact solution must have maximum objective
    expect_equal(which.max(obj_values), 1)
  }
})
