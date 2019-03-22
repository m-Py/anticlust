
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
    anticlusters_exact <- anticlustering(features, n_clusters, method = "exact", standardize = FALSE,
                                         preclustering = FALSE)
    anticlusters_heuristic <- anticlustering(features, n_clusters,
                                             method = "sampling",
                                             standardize = FALSE,
                                             nrep = 5)
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


test_that("all argument combinations run through", {
  conditions <- expand.grid(preclustering = c(TRUE, FALSE),
                            method = c("exact", "annealing", "sampling"))
  # Set up matrix to store the objective values obtained by different methods
  storage <- matrix(ncol = 3, nrow = 2)
  colnames(storage) <- c("exact", "annealing", "sampling")
  rownames(storage) <- c("preclustering", "no_preclustering")

  criterion <- "distance"
  n_elements <- 12
  features <- matrix(round(rnorm(n_elements * 2)), ncol = 2)
  n_anticlusters <- 2

  for (i in 1:nrow(conditions)) {
    method <- conditions$method[i]
    preclustering <- conditions$preclustering[i]
    anticlusters <- anticlustering(features, n_anticlusters, criterion,
                                   method, preclustering, standardize = FALSE)
    obj <- get_objective(features, anticlusters, objective = criterion)
    rowname <- ifelse(preclustering, "preclustering", "no_preclustering")
    storage[rowname, method] <- obj
  }
  ## Exact solution must be best:
  expect_equal(all(round(storage["no_preclustering", "exact"], 10) >= round(c(storage), 10)), TRUE)
})
