
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
                              method = "ilp", standardize = FALSE)
    ac_dist <- anticlustering(distances = distances, K = K,
                              preclustering = FALSE,
                              method = "ilp", standardize = FALSE)
    expect_equal(distance_objective_(ac_feat, dist(features)),
                 distance_objective_(ac_feat, distances))
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
                              method = "ilp", standardize = FALSE)
    ac_dist <- anticlustering(distances = distances, K = K,
                              preclustering = TRUE,
                              method = "ilp", standardize = FALSE)
    expect_equal(distance_objective_(ac_feat, dist(features)),
                 distance_objective_(ac_feat, distances))
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

    expect_equal(distance_objective_(ac_feat, dist(features)),
                 distance_objective_(ac_feat, distances))
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

    ## Use a fixed seed to compare the random sampling method based on
    ## features and distance input
    rnd_seed <- sample(10000, size = 1)

    set.seed(rnd_seed)
    ac_feat <- anticlustering(features, K = K, preclustering = FALSE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)
    set.seed(rnd_seed)
    ac_dist <- anticlustering(distances = distances, K = K, preclustering = FALSE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)

    expect_equal(distance_objective_(ac_feat, dist(features)),
                 distance_objective_(ac_feat, distances))
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

    ## Use a fixed seed to compare the random sampling method based on
    ## features and distance input
    rnd_seed <- sample(10000, size = 1)

    set.seed(rnd_seed)
    ## does not work with distance criterion
    ac_feat <- anticlustering(features, K = K, preclustering = TRUE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)
    set.seed(rnd_seed)
    ac_dist <- anticlustering(distances = distances, K = K, preclustering = TRUE,
                              method = "heuristic", standardize = FALSE,
                              nrep = 100)

    expect_equal(distance_objective_(ac_feat, dist(features)),
                 distance_objective_(ac_feat, distances))
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

    expect_equal(distance_objective_(ac_feat, dist(features)),
                 distance_objective_(ac_feat, distances))
  }
})
