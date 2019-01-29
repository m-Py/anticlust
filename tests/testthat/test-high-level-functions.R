
context("Test high-level functions for equal sized clustering and anticlustering")
library("anticlust")

test_that("high level equal sized clustering function runs through", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    n_clusters <- conditions[k, "p"]
    n_elements <- n_clusters * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    clusters_exact <- clustering(features, n_clusters, method = "exact", standardize = FALSE)
    clusters_heuristic <- clustering(features, n_clusters, method = "heuristic", standardize = FALSE)
    ## Check that output is valid
    expect_equal(legal_number_of_clusters(features, clusters_exact), NULL)
    expect_equal(legal_number_of_clusters(features, clusters_heuristic), NULL)

    ## Assert that exact solution has lowest objective (for distance
    ## criterion), allowing for numeric imprecision of ILP solver
    obj_exact     <- get_objective(features, clusters_exact, "distance")
    obj_heuristic <- get_objective(features, clusters_heuristic, "distance")
    expect_equal(round(obj_exact, 10) <= round(obj_heuristic, 10), TRUE)
  }
})



test_that("high level anticlustering function runs through", {
  conditions <- expand.grid(m = 1:4, p = 2:3)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    n_clusters <- conditions[k, "p"]
    n_elements <- n_clusters * 3 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    anticlusters_exact <- anticlustering(features, n_clusters, method = "exact", standardize = FALSE)
    anticlusters_heuristic <- anticlustering(features, n_clusters,
                                             method = "random", standardize = FALSE)
    ## Check that output is valid
    expect_equal(legal_number_of_clusters(features, anticlusters_exact), NULL)
    expect_equal(legal_number_of_clusters(features, anticlusters_heuristic), NULL)
    ## Assert that exact solution has highest objective (for distance
    ## criterion), allowing for numeric imprecision of ILP solver
    obj_exact     <- get_objective(features, anticlusters_exact, "distance")
    obj_heuristic <- get_objective(features, anticlusters_heuristic, "distance")
    expect_equal(round(obj_exact, 10) >= round(obj_heuristic, 10), TRUE)
  }
})
