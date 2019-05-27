
context("Test heuristic methods for (anti)clustering")
library("anticlust")

test_that("heuristic clustering works for different inputs", {
  for (m in 1:4) {
    m_features <- m
    ## vary number of anticlusters
    for (p in 2:5) {
      p_anticlusters <- p
      n_elements <- p * 3 # n-elements must be multiplier of p
      features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
      clusters <- equal_sized_kmeans(features, p)
      expect_equal(legal_number_of_clusters(features, clusters), NULL)
    }
  }
})

test_that("heuristic clustering produces expected output", {
  ## Vary number of features
  for (m in 1:4) {
    m_features <- m
    ## Vary number of clusters
    for (p in 2:4) {
      p_anticlusters <- p
      n_elements <- p * 5 # n-elements must be multiplier of p
      features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
      cluster_init <- kmeans(features, p_anticlusters)$cl
      centers  <- cluster_centers(features, cluster_init)
      clust_dist <- dist_from_centers(features, centers, squared = TRUE)
      clusters <- heuristic_cluster_assignment(clust_dist)
      ## Assert that each element is assigned only once
      expect_equal(sum(duplicated(clusters)), 0)
      ## Check dimensionality of output
      expect_equal(ncol(clusters), p)
      expect_equal(nrow(clusters), n_elements / p)

      ## Test the cluster_to_long method (is the correct order preserved?)
      long_clusters <- clusters_to_long(clusters)
      expect_equal(legal_number_of_clusters(features, long_clusters), NULL)
      for (i in 1:nrow(clusters)) {
        for (j in 1:ncol(clusters)) {
          index <- clusters[i, j]
          ## The element at this index in long_cluster must correspond
          ## to the column in matrix `cluster` (both encode the cluster)
          expect_equal(long_clusters[index], j)
        }
      }
    }
  }
})

test_that("heuristic anticlustering produces expected output", {
  conditions <- expand.grid(m = 1:4, p = 2:4, objective = c("distance", "variance"),
                            preclustering = c(TRUE, FALSE))
  for (i in 1:nrow(conditions)) {
    m_features <- conditions[i, "m"]
    p_anticlusters <- conditions[i, "p"]
    n_elements <- p_anticlusters * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)

    ## First: heuristic preclustering
    preclusters <- NULL
    if (conditions$preclustering[i] == TRUE) {
      n_preclusters <- n_elements / p_anticlusters
      preclusters <- equal_sized_kmeans(features, n_preclusters)
      ## Legal number of preclusters?
      expect_equal(legal_number_of_clusters(features, preclusters), NULL)
      ## Expected number of preclusters?
      expect_equal(as.numeric(table(preclusters)[1]), n_elements / n_preclusters)
    }

    ## Now anticlustering:
    # Test parallel and non-parallel anticlustering:
    for (parallel in c(TRUE, FALSE)) {
      anticlusters <- heuristic_anticlustering(as.matrix(features), p_anticlusters,
                                               preclusters, nrep = 100,
                                               objective = conditions$objective[i],
                                               categories = NULL, parallelize = parallel,
                                               seed = NULL)
      ## Legal number of anticlusters?
      expect_equal(legal_number_of_clusters(features, anticlusters), NULL)
      ## Expected number of anticlusters?
      expect_equal(as.numeric(table(anticlusters)[1]), n_elements / p_anticlusters)
    }
  }
})
