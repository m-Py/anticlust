
# Wrap the balanced clustering function to allow for some imbalancing
# @param features
# @param K the Number of anticlusters (! -> the number of clusters returned is ~ nrow(features) / K)

imbalanced_preclustering <- function(features, K) {
  features <- as.matrix(features)
  N <- nrow(features)
  preclusters <- rep(NA, N)
  # only select as many data as can be clustered into balanced clusters
  subsetted <- features[1:(N - (N %% K)), , drop = FALSE]
  preclusters[1:nrow(subsetted)] <- balanced_clustering(
    features,
    K = nrow(subsetted) / K
  )
  # full some clusters at random
  if (sum(is.na(preclusters)) > 0) {
    preclusters[is.na(preclusters)] <- sample(
      1:max(preclusters, na.rm = TRUE), 
      size = sum(is.na(preclusters))
    )
  }
  preclusters
}
