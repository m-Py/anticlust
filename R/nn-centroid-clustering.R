
# nearest neighbor centroid clustering

# data = feature data frame
# K = size of the small groups
nn_centroid_clustering <- function(data, K) {
  data <- as.matrix(data)
  # compute the distances between all elements and overall cluster center
  centroid <- t(as.matrix(colMeans(data)))
  distances <- c(dist_from_centers(data, centroid, squared = FALSE))
  counter <- 1
  idx <- 1:nrow(data)
  clusters <- rep(NA, nrow(data))
  
  while (nrow(data) > 0) {
    # get element with maximum distance to centroid
    max_away <- which.max(distances)[1]
    # compute nearest neighbors for element that is furthest away
    clustered <- nn2(data, data[max_away, , drop = FALSE], K)$nn.idx
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
