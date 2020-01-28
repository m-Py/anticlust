
# nearest neighbor centroid clustering

# data = feature data frame
# K = size of the small groups (not the number of clusters!)
# groups = vector of length nrow(data); if passed, K-partite matching is conducted,
# i.e., each element is clustered with elements that are in *other* groups than 
# the element itself. Consists of integers 1, ..., K
nn_centroid_clustering <- function(
  data, K, 
  groups = NULL, 
  match_extreme_first = TRUE, 
  target_group = FALSE
) {
  data <- as.matrix(data)
  N <- nrow(data)
  distances <- distances_from_centroid(data)

  if (argument_exists(groups)) {
    K <- length(unique(groups))
  }

  # some book-keeping variables for the loop
  cluster_id <- 1
  idx <- 1:nrow(data)
  clusters <- rep(NA, nrow(data))
  while (data_fits(data, K) && no_group_empty(groups, K)) {
    # Get target item
    target <- get_target(distances, groups, target_group, match_extreme_first)
    # Compute nearest neighbors for target element
    clustered <- get_nearest_neighbours(data, target, K, groups)
    clusters[idx[clustered]] <- cluster_id
    # Remove matched elements:
    data <- subset_data_matrix(data, -clustered)
    distances <- distances[-clustered]
    groups <- groups[-clustered]
    idx  <- idx[-clustered]
    # increase cluster id for next elements
    cluster_id <- cluster_id + 1
  }
  order_cluster_vector(clusters)
}

# simple test: does a data frame still fit into pieces of K units
data_fits <- function(data, K)  {
  nrow(data) > (K-1)
}

# for k-partite matching: test if no group is empty 
no_group_empty <- function(groups, K) {
  !(argument_exists(groups) && any(!(1:K %in% groups)))
}

# Find target element for which neighbours are sought
# param distances: distances to a cluster centroid
# param groups: grouping vector (%in% 1, ..., K)
# param target_group: id of the group from which target is selected 
#    (may be FALSE instead of an ID)
# param match_extreme_first: boolean
get_target <- function(distances, groups, target_group, match_extreme_first) {
  if (!target_group) {
    # no target group specified: select an item from all elements
    if (match_extreme_first) {
      return(which.max(distances))
    } else {
      return(which.min(distances))
    }
  }
  # return a member from the target group as target
  ids_target <- which(groups == target_group)
  ordered_distances <- order(distances, decreasing = match_extreme_first)
  ordered_distances[ordered_distances %in% ids_target][1]
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

# Subset a distance or feature data matrix (not knowing which one)
subset_data_matrix <- function(data, selection) {
  if (is_distance_matrix(data)) {
    data <- data[selection, selection, drop = FALSE]
  } else {
    data <- data[selection, , drop = FALSE]
  }
  data
}
