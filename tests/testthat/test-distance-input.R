
context("Test distance input")
library("anticlust")

test_that("distance input works for exact ILP", {
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)

    ac_feat <- anticlustering(features, K = K, preclustering = FALSE,
                              method = "exact", standardize = FALSE)
    ac_dist <- anticlustering(distances = distances, K = K,
                              preclustering = FALSE,
                              method = "exact", standardize = FALSE)
    expect_equal(obj_value_distance(features, ac_feat),
                 distance_objective(distances, ac_feat))
    expect_equal(obj_value_distance(features, ac_dist),
                 distance_objective(distances, ac_dist))
  }
})

test_that("distance input works for precluster ILP", {
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)

    ac_feat <- anticlustering(features, K = K, preclustering = TRUE,
                              method = "exact", standardize = FALSE)
    ac_dist <- anticlustering(distances = distances, K = K,
                              preclustering = TRUE,
                              method = "exact", standardize = FALSE)
    expect_equal(obj_value_distance(features, ac_feat),
                 distance_objective(distances, ac_feat))
    expect_equal(obj_value_distance(features, ac_dist),
                 distance_objective(distances, ac_dist))
  }
})

test_that("distance input works for complete enumeration", {
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)

    ac_feat <- enum_anticlustering(features, K = K)
    ac_dist <- enum_anticlustering(distances = distances, K = K)

    expect_equal(obj_value_distance(features, ac_feat),
                 distance_objective(distances, ac_feat))
    expect_equal(obj_value_distance(features, ac_dist),
                 distance_objective(distances, ac_dist))
  }
})

test_that("distance input works for heuristic without preclustering", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)

    ac_feat <- anticlustering(features, K = K, preclustering = FALSE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)
    ac_dist <- anticlustering(distances = distances, K = K, preclustering = FALSE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)

    expect_equal(obj_value_distance(features, ac_feat),
                 distance_objective(distances, ac_feat))
    expect_equal(obj_value_distance(features, ac_dist),
                 distance_objective(distances, ac_dist))
  }
})

test_that("distance input works for heuristic with preclustering", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)

    ## does not work with distance criterion
    ac_feat <- anticlustering(features, K = K, preclustering = TRUE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)
    ac_dist <- anticlustering(distances = distances, K = K, preclustering = TRUE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)

    expect_equal(obj_value_distance(features, ac_feat),
                 distance_objective(distances, ac_feat))
    expect_equal(obj_value_distance(features, ac_dist),
                 distance_objective(distances, ac_dist))
  }
})


test_that("distance input works for clustering function, heuristic method", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- dist(features)

    ## does not work with distance criterion
    ac_feat <- balanced_clustering(features, K = K, method = "heuristic", standardize = FALSE)
    ac_dist <- balanced_clustering(distances = distances, K = K, method = "heuristic", standardize = FALSE)

    expect_equal(obj_value_distance(features, ac_feat),
                 distance_objective(distances, ac_feat))
    expect_equal(obj_value_distance(features, ac_dist),
                 distance_objective(distances, ac_dist))
  }
})
