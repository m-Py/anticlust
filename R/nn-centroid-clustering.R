
# nearest neighbor centroid clustering

# data = feature data frame
# K = size of the small groups
nn_centroid_clustering <- function(data, K) {
  data <- as.matrix(data)
  if (sum(is.na(data)) != 0) {
    data <- impute_mean_na(data)
  }
  
  distances <- distances_from_centroid(data)
  counter <- 1
  idx <- 1:nrow(data)
  clusters <- rep(NA, nrow(data))
  
  while (nrow(data) > 0) {
    # compute nearest neighbors for element that is furthest away
    clustered <- get_nearest_neighbours(data, which.max(distances), K)
    clusters[idx[clustered]] <- counter
    data <- data[-clustered, , drop = FALSE]
    distances <- distances[-clustered]
    idx  <- idx[-clustered]
    counter <- counter + 1
  }
  if (sum(is.na(clusters)) > 0) {
    stop("something is wrong")
  }
  clusters
}

# Get nearest neighbours for a current element
# param data: the data
# param may_away: the index of the element that is furthest away from 
#   the cluster centroid and for which nearest neighbors are sought
# param K: The number of clusters / Nearest neighbours 
# return: The indices of the current element as vector. !! The first
#   index is the element itself !!
get_nearest_neighbours <- function(data, max_away, K) {
  nn2(data, data[max_away, , drop = FALSE], K)$nn.idx
}

# Compute the distances from centroid of a data set
distances_from_centroid <- function(data) {
  centroid <- t(as.matrix(colMeans(data)))
  c(dist_from_centers(data, centroid, squared = FALSE))
}

# for NA values use, mean imputation
impute_mean_na <- function(data) {
  warning("I used mean imputation for NA values; you really should ",
          "consider if you are okay with having missing ",
          "values in your data!")
  apply(data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
}
