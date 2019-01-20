library("anticlust")

context("KMeans-Maximizing")

## What data structures does k-means use?

# 1. Data matrix of features - row is element, column is feature.
#    n rows, m features
# 2. A vector of anticluster assignments
#    n elements, p unique elements (= number of different anticlusters)
# 3. Matrix of anticluster centers. Each row corresponds to a
#    cluster and each column corresponds to a feature.
#    p rows, m columns
# 4. Matrix, distances from points to anticluster centers,
#    p columns, n rows

test_that("computation of cluster centers has correct output", {
  ## Vary number of features
  for (m in 1:4) {
    m_features <- m
    ## vary number of anticlusters
    for (p in 2:5) {
      p_anticlusters <- p
      n_elements <- p * 3 # n-elements must be multiplier of p
      features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
      anticlusters <- rep(1:p_anticlusters, n_elements/p_anticlusters)
      centers  <- cluster_centers(features, anticlusters)
      ## Chech output of cluster_centers
      expect_equal(class(centers), "matrix")
      # Matrix of distances should have p rows
      expect_equal(dim(centers)[1], p_anticlusters)
      # Matrix of distances should have m columns
      expect_equal(dim(centers)[2], m_features)
      ## Check that means are computed correctly for all features
      for (i in 1:m) {
        expect_equal(all(centers[, i] == tapply(features[, i], anticlusters, mean)), TRUE)
      }
    }
  }
})

## TODO next: implement test for computation of distances from cluster centers
