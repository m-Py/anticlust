
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
# apart from any other data points
remove_outliers <- function(data, equalize, K) {
  N <- nrow(data)
  if (N %% K == 0) {
    return(data)
  }
  distances <- as.matrix(dist(data[, equalize]))
  diag(distances) <- Inf
  minima <- apply(distances, 1, min)
  minima <- cbind(minima, 1:N)
  minima <- sort_by_col(minima, 1)
  # discard elements whose nearest neighbor is farthest away
  minima <- minima[1:(N - (N %% K)), , drop = FALSE]
  # obtain original order
  minima <- sort_by_col(minima, 2)
  data[minima[, 2], , drop = FALSE]
}
