
library("anticlust")

# heuristic clustering works for different inputs
for (m in 1:4) {
  m_features <- m
  ## vary number of anticlusters
  for (p in 2:5) {
    n_elements <- p * 3 # n-elements must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    clusters <- balanced_clustering(features, K = p)
    expect_equal(anticlust:::legal_number_of_clusters(features, clusters), NULL)
    clusters <- balanced_clustering(dist(features), K = p)
    expect_equal(anticlust:::legal_number_of_clusters(features, clusters), NULL)
  }
}

# heuristic anticlustering produces expected output
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
    preclusters <- balanced_clustering(features, K = n_preclusters)
    ## Legal number of preclusters?
    expect_equal(anticlust:::legal_number_of_clusters(features, preclusters), NULL)
    ## Expected number of preclusters?
    expect_equal(as.numeric(table(preclusters)[1]), n_elements / n_preclusters)
  }
  
  ## Now anticlustering:
  obj_function <- ifelse(conditions$objective[i] == "distance",
                         diversity_objective, variance_objective)
  anticlusters <- anticlustering(
    features,
    K = p_anticlusters,
    categories = preclusters,
    objective = obj_function
  )
  ## Legal number of anticlusters?
  expect_equal(anticlust:::legal_number_of_clusters(features, anticlusters), NULL)
  ## Expected number of anticlusters?
  expect_equal(as.numeric(table(anticlusters)[1]), n_elements / p_anticlusters)
  if (conditions$preclustering[i] == TRUE) {
    # Are no preclustered elements together?
    tab <- table(anticlusters, preclusters)
    expect_equal(all(tab == 1), TRUE)
  }
}
