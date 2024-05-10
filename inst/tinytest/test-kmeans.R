library("anticlust")


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


## A standard (or squared) euclidean distance between two data points
## Currently, this function is only used in the test files.
euc_dist <- function(x1, x2, squared = FALSE) {
  if (squared)
    return(sum((x1 - x2)^2))
  sqrt(sum((x1 - x2)^2))
}


# computation of cluster centers and distances to centers is correct
## Vary number of features
for (m in 1:4) {
  m_features <- m
  ## vary number of anticlusters
  for (p in 2:5) {
    p_anticlusters <- p
    n_elements <- p * 3 # n-elements must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    anticlusters <- rep(1:p_anticlusters, n_elements/p_anticlusters)
    centers  <- anticlust:::cluster_centers(features, anticlusters)
    ## Check output of cluster_centers
    expect_equal(inherits(centers, "matrix"), TRUE)
    # Matrix of distances should have p rows
    expect_equal(dim(centers)[1], p_anticlusters)
    # Matrix of distances should have m columns
    expect_equal(dim(centers)[2], m_features)
    ## Check that means are computed correctly for all features
    for (i in 1:m) {
      expect_equal(as.numeric(centers[, i]),
                   as.numeric(tapply(features[, i], anticlusters, mean)))
    }
    
    ## Now also check if distances to cluster centers are computed correctly
    distances <- anticlust:::dist_from_centers(features, centers, squared = FALSE)
    expect_equal(inherits(centers, "matrix"), TRUE)
    expect_equal(dim(distances)[1], n_elements)
    expect_equal(dim(distances)[2], p_anticlusters)
    ## Check distance output for all points
    for (f in 1:n_elements) {
      for (c in 1:p_anticlusters) {
        expect_equal(distances[f, c], euc_dist(features[f, ], centers[c, ]))
      }
    }
  }
}
