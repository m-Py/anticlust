
# nearest neighbor centroid clustering

# data = feature data frame
# K = size of the small groups (not the number of clusters!)
# groups = vector of length nrow(data); if passed, K-partite matching is conducted,
# i.e., each element is clustered with elements that are in *other* groups than 
# the element itself
nn_centroid_clustering <- function(data, K, groups = NULL, match_extreme_first = TRUE) {
  data <- as.matrix(data)
  N <- nrow(data)
  distances <- distances_from_centroid(data)

  if (argument_exists(groups)) {
    # smallest group determines what matches are possible
    smallest_group <- which.min(table(groups))
    K <- length(unique(groups))
  }
  
  # some book-keeping variables for the loop
  counter <- 1
  idx <- 1:nrow(data)
  clusters <- rep(NA, nrow(data))
  
  while (nrow(data) > (K-1)) {
    # compute nearest neighbors for a target element
    target <- get_target(distances, groups, smallest_group, match_extreme_first)
    clustered <- get_nearest_neighbours(data, target, K, groups)
    clusters[idx[clustered]] <- counter
    if (is_distance_matrix(data)) {
      data <- data[-clustered, -clustered, drop = FALSE]
    } else {
      data <- data[-clustered, , drop = FALSE]
    }
    distances <- distances[-clustered]
    groups <- groups[-clustered]
    idx  <- idx[-clustered]
    counter <- counter + 1
    # for bipartite matching: stop as soon the smallest group is matched
    # with elements from the other groups
    if (argument_exists(groups)) {
      if (sum(groups == smallest_group) == 0) {
        break
      }
    }
  }
  # for subset matching: fill the remaining clusters: 
  order_cluster_vector(clusters)
}

# Find target element for which neighbours are sought
# param distances: distances to a cluster centroid
# param groups: groupingg vector
# param smallest_group: id of the group having the fewest members
get_target <- function(distances, groups, smallest_group, match_extreme_first) {
  # if bipartite subset selection is required: 
  # select a member from the smallest group
  if (argument_exists(groups)) {
    group_sizes <- table(groups)
    if (any(group_sizes != group_sizes[1])) {
      # return the most extreme member from the smallest group as target
      ids_smallest <- which(groups == smallest_group)
      ordered_distances <- order(distances, decreasing = match_extreme_first)
      return(ordered_distances[ordered_distances %in% ids_smallest][1])
    }
  }
  # otherwise: select an item from all elements
  if (match_extreme_first) {
    return(which.max(distances))
  }
  which.min(distances)
}


# Get nearest neighbours for a current element
# param data: the data
# param target: the index of the element for which nearest neighbors are sought
# param K: The number of nearest neighbours (including the element itself!)
# param groups = vector of length nrow(data); elements are clustered with elements in *other* groups 
# return: The indices of the current element as vector. !! The first
#   index is the element itself !!
get_nearest_neighbours <- function(data, target, K, groups) {
  
  # k-partite clustering
  if (argument_exists(groups)) {
    groups <- to_numeric(groups)
    current_group <- groups[target]
    other_groups <- (1:K)[-current_group]
    current_element <- data[target, , drop = FALSE]
    # get nearest neighbour in each other group
    nns <- rep(NA, K-1)
    for (i in seq_along(other_groups)) {
      if (is_distance_matrix(data)) {
        # get indices of elements in other groups
        idx_other <- which(groups == other_groups[i])
        most_similar <- order(current_element)
        nns[i] <- most_similar[most_similar %in% idx_other][1]
      } else {
        tmp_data <- data[groups == other_groups[i], , drop = FALSE]
        nn <- c(nn2(tmp_data, current_element, 1)$nn.idx)
        nns[i] <- which(groups == other_groups[i])[nn]
      }
    }
    nns <- c(target, nns)
    return(nns)
  }
  
  # normal clustering - no grouping restrictions
  if (is_distance_matrix(data)) {
    nns <- order(data[target, ])[1:K]
  } else {
    nns <- nn2(data, data[target, , drop = FALSE], K)$nn.idx
  }
  nns
}

# Compute the distances from centroid of a data set.
# Centroid is either computed in euclidean space or taken as a central element
distances_from_centroid <- function(x) {
  if (is_distance_matrix(x)) {
    distances <- as.matrix(x)
    maxima <- apply(distances, 1, max)
    centroid <- which.min(maxima)
    return(distances[centroid, ])
  }
  centroid <- t(as.matrix(colMeans(x)))
  c(dist_from_centers(x, centroid, squared = FALSE))
}


#' Is a matrix a legal distance matrix
#'
#' @param m a Matrix
#' @return TRUE if `m` is distance matrix, FALSE otherwise
#' @noRd
is_distance_matrix <- function(m) {
  m <- as.matrix(m)
  if (nrow(m) != ncol(m)) {
    return(FALSE)
  }
  lower <- m[lower.tri(m)]
  m <- t(m)
  upper <- m[lower.tri(m)]
  all(lower == upper)
}
