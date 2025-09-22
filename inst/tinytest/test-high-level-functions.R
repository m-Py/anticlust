
library("anticlust")

# high level equal sized clustering function runs through
conditions <- expand.grid(m = 1:4, p = 2:4)
for (k in 1:nrow(conditions)) {
  m_features <- conditions[k, "m"]
  n_clusters <- conditions[k, "p"]
  n_elements <- n_clusters * 5 # n must be multiplier of p
  features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
  clusters_exact <- balanced_clustering(features, K = n_clusters, method = "ilp")
  clusters_heuristic <- balanced_clustering(features, K = n_clusters,
                                            method = "centroid")
  ## Check that output is valid
  expect_equal(anticlust:::legal_number_of_clusters(features, clusters_exact), NULL)
  expect_equal(anticlust:::legal_number_of_clusters(features, clusters_heuristic), NULL)
  
  ## Assert that exact solution has lowest objective (for distance
  ## criterion), allowing for numeric imprecision of ILP solver
  obj_exact     <- anticlust:::diversity_objective_(clusters_exact, features)
  obj_heuristic <- anticlust:::diversity_objective_(clusters_heuristic, features)
  expect_equal(round(obj_exact, 10) <= round(obj_heuristic, 10), TRUE)
}


# high level anticlustering function runs through
conditions <- expand.grid(m = 1:4, p = 2:3)
for (k in 1:nrow(conditions)) {
  m_features <- conditions[k, "m"]
  n_clusters <- conditions[k, "p"]
  n_elements <- n_clusters * 3 # n must be multiplier of p
  features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
  anticlusters_exact <- anticlustering(features, K = n_clusters,
                                       method = "ilp",
                                       preclustering = FALSE)
  anticlusters_heuristic <- anticlustering(features,
                                           K = n_clusters)
  ## Check that output is valid
  expect_equal(anticlust:::legal_number_of_clusters(features, anticlusters_exact), NULL)
  expect_equal(anticlust:::legal_number_of_clusters(features, anticlusters_heuristic), NULL)
  ## Assert that exact solution has highest objective (for distance
  ## criterion), allowing for numeric imprecision of ILP solver
  obj_exact     <- anticlust:::diversity_objective_(anticlusters_exact, features)
  obj_heuristic <- anticlust:::diversity_objective_(anticlusters_heuristic, features)
  expect_equal(round(obj_exact, 10) >= round(obj_heuristic, 10), TRUE)
}


# all argument combinations run through
conditions <- expand.grid(preclustering = c(TRUE, FALSE),
                          method = c("ilp", "exchange"))
# Set up matrix to store the objective values obtained by different methods
storage <- matrix(ncol = 2, nrow = 2)
colnames(storage) <- c("ilp", "exchange")
rownames(storage) <- c("preclustering", "no_preclustering")

criterion <- "distance"
n_elements <- 12
features <- matrix(rnorm(n_elements * 2), ncol = 2)
n_anticlusters <- 2

for (i in 1:nrow(conditions)) {
  method <- conditions$method[i]
  preclustering <- conditions$preclustering[i]
  anticlusters <- anticlustering(features, K = n_anticlusters,
                                 objective = criterion,
                                 method = method,
                                 preclustering = preclustering)
  obj <- anticlust:::diversity_objective_(anticlusters, features)
  rowname <- ifelse(preclustering, "preclustering", "no_preclustering")
  storage[rowname, method] <- obj
}
## Exact solution must be best:
expect_equal(all(round(storage["no_preclustering", "ilp"], 10) >= round(c(storage), 10)), TRUE)

