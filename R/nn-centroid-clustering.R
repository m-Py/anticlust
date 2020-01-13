
# nearest neighbor centroid clustering

# data = feature data frame
# K = size of the small groups
# groups = vector of length nrow(data); if passed, K-partite matching is conducted,
# i.e., each element is clustered with elements that are in *other* groups than 
# the element itself
nn_centroid_clustering <- function(data, K, groups = NULL) {
  data <- as.matrix(data)

  distances <- distances_from_centroid(data)
  counter <- 1
  idx <- 1:nrow(data)
  clusters <- rep(NA, nrow(data))
  
  while (nrow(data) > 0) {
    # compute nearest neighbors for element that is furthest away
    max_away <- which.max(distances)
    clustered <- get_nearest_neighbours(data, max_away, K, groups)
    # as soon as no more neighbours are found: exit
    if (any(is.na(clustered))) {
      break
    }
    clusters[idx[clustered]] <- counter
    data <- data[-clustered, , drop = FALSE]
    distances <- distances[-clustered]
    groups <- groups[-clustered]
    idx  <- idx[-clustered]
    counter <- counter + 1
  }
  clusters
}

# Get nearest neighbours for a current element
# param data: the data
# param may_away: the index of the element that is furthest away from 
#   the cluster centroid and for which nearest neighbors are sought
# param K: The number of clusters / nearest neighbours
# param groups = vector of length nrow(data); elements are clustered with elements in *other* groups 
# return: The indices of the current element as vector. !! The first
#   index is the element itself !!
get_nearest_neighbours <- function(data, max_away, K, groups) {
  if (argument_exists(groups)) {
    groups <- to_numeric(groups)
    current_group <- groups[max_away]
    other_groups <- (1:K)[-current_group]
    current_element <- data[max_away, , drop = FALSE]
    
    # get nearest neighbour in each other group
    nns <- rep(NA, K-1)
    for (i in seq_along(other_groups)) {
      tmp_data <- data[groups == other_groups[i], , drop = FALSE]
      
      ## check if a neighbor can be found at all; otherwise return NA!
      if (nrow(tmp_data) > 0) {
        nn <- c(nn2(tmp_data, current_element, 1)$nn.idx)
        nns[i] <- which(groups == other_groups[i])[nn]
      } else {
        nns[i] <- NA
      }
    }
    nns <- c(max_away, nns)
    return(nns)
  }
  nn2(data, data[max_away, , drop = FALSE], K)$nn.idx
}

# Compute the distances from centroid of a data set
distances_from_centroid <- function(data) {
  centroid <- t(as.matrix(colMeans(data)))
  c(dist_from_centers(data, centroid, squared = FALSE))
}

}
