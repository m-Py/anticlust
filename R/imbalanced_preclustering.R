
# Wrap the balanced clustering function to allow for some imbalancing
# @param features
# @param K the Number of anticlusters (! -> the number of clusters returned is ~ nrow(features) / K)

imbalanced_preclustering <- function(features, K) {
  features <- as.matrix(features)
  N <- nrow(features)
  preclusters <- rep(NA, N)
  # only select as many data as can be clustered into balanced clusters
  subsetted <- remove_outliers(features, 1:ncol(features), K)
  preclusters[1:nrow(subsetted)] <- balanced_clustering(
    subsetted,
    K = nrow(subsetted) / K
  )
  # the remainders get into the last cluster
  if (sum(is.na(preclusters)) > 0) {
    preclusters[is.na(preclusters)] <- max(preclusters, na.rm = TRUE) + 1
  }
  preclusters
}
