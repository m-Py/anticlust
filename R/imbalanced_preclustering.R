
# Wrap the balanced clustering function to allow for some imbalancing
# @param features
# @param K the Number of anticlusters (! -> the number of clusters returned is ~ nrow(features) / K)

imbalanced_preclustering <- function(features, K) {
  features <- data.frame(features)
  N <- nrow(features)
  features$preclusters_id <- 1:N
  preclusters <- rep(NA, N)
  # only select as many data as can be clustered into balanced clusters
  subsetted <- remove_outliers(features, K)
  # test if only one cluster can be filled:
  if (nrow(subsetted) == K)  {
    preclusters[subsetted$preclusters_id] <- 1
  } else {
    cl <- nn_centroid_clustering(
      subsetted[, -ncol(subsetted)],
      K = K
    )
    preclusters[subsetted$preclusters_id] <- cl
  }
  # the remainders get into the last cluster
  if (sum(is.na(preclusters)) > 0) {
    preclusters[is.na(preclusters)] <- max(preclusters, na.rm = TRUE) + 1
  }
  preclusters
}


# Function to remove elements such that the remaining elements can fit 
# into clusters of size K. Removes the data points that are furthest 
# apart from the data centroid
remove_outliers <- function(data, K) {
  N <- nrow(data)
  if (N %% K == 0) {
    return(data)
  }
  # define outliers as being furthest away from data centroid
  distances <- distances_from_centroid(data)
  # get indices of elements that are no outliers
  idx <- order(distances)[1:(N - (N %% K))]
  # return these elements
  data[idx, , drop = FALSE]
}
