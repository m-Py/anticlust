
context("Test distance input")
library("anticlust")

test_that("distance input works for exact ILP", {
  skip_on_cran()
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- as.matrix(dist(features))

    ac_feat <- anticlustering(features, K = K, preclustering = FALSE,
                              method = "ilp")
    ac_dist <- anticlustering(distances, K = K,
                              preclustering = FALSE,
                              method = "ilp")
    expect_equal(diversity_objective_(ac_feat, distances),
                 diversity_objective_(ac_feat, features))
    expect_equal(diversity_objective_(ac_dist, distances),
                 diversity_objective_(ac_dist, features))
    expect_equal(diversity_objective_(ac_feat, distances),
                 diversity_objective_(ac_dist, distances))
  }
})

test_that("distance input works for precluster ILP", {
  skip_on_cran()
  conditions <- expand.grid(m = 1:4, p = 2)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- as.matrix(dist(features))

    ac_feat <- anticlustering(features, K = K, preclustering = TRUE,
                              method = "ilp")
    ac_dist <- anticlustering(distances, K = K,
                              preclustering = TRUE,
                              method = "ilp")
    expect_equal(diversity_objective_(ac_feat, distances),
                 diversity_objective_(ac_feat, features))
    expect_equal(diversity_objective_(ac_dist, distances),
                 diversity_objective_(ac_dist, features))
    expect_equal(diversity_objective_(ac_feat, distances),
                 diversity_objective_(ac_dist, distances))
  }
})

test_that("distance input works for exchange method", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    K <- conditions[k, "p"]
    n_elements <- K * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- as.matrix(dist(features))

    ## Use a fixed seed to compare the random sampling method based on
    ## features and distance input
    rnd_seed <- sample(10000, size = 1)

    set.seed(rnd_seed)
    ac_feat <- anticlustering(features, K = K)
    set.seed(rnd_seed)
    ac_dist <- anticlustering(distances, K = K)

    expect_equal(diversity_objective_(ac_feat, distances),
                 diversity_objective_(ac_feat, features))
    expect_equal(diversity_objective_(ac_dist, distances),
                 diversity_objective_(ac_dist, features))
    expect_equal(diversity_objective_(ac_feat, distances),
                 diversity_objective_(ac_dist, distances))
  }
})
