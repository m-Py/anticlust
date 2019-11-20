
# Wrap the balanced clustering function to allow for some imbalancing
# @param features
# @param K the Number of anticlusters (! -> the number of clusters returned is ~ nrow(features) / K)

imbalanced_preclustering <- function(features, K) {
  features <- data.frame(features)
  N <- nrow(features)
  features$preclusters_id <- 1:N
  preclusters <- rep(NA, N)
  # only select as many data as can be clustered into balanced clusters
  subsetted <- remove_outliers(features, 1:ncol(features), K)
  # test if only one cluster can be filled:
  if (nrow(subsetted) == K)  {
    preclusters[subsetted$preclusters_id] <- 1
  } else {
    preclusters[subsetted$preclusters_id] <- balanced_clustering(
      subsetted,
      K = nrow(subsetted) / K
    )
  }
  # the remainders get into the last cluster
  if (sum(is.na(preclusters)) > 0) {
    preclusters[is.na(preclusters)] <- max(preclusters, na.rm = TRUE) + 1
  }
  preclusters
}
