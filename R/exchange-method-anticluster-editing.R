
#' Solve anticluster editing using a fast exchange method
#'
#' @param data An N x N dissimilarity matrix or N x M features matrix.
#' @param K The number of cluster or an initial cluster assignment
#' @param categories A vector representing categorical constraints
#' @param preclusters A vector representing preclustering constraints
#'
#' @return The anticluster assignment
#'
#' @noRd
#'

fast_exchange_dist <- function(data, K, categories) {
  distances <- convert_to_distances(data)  
  clusters <- initialize_clusters(nrow(distances), K, categories)
  best_total <- diversity_objective_(clusters, distances)
  distances[upper.tri(distances)] <- 0
  diag(distances) <- 0
  ## Matrix for indexing the elements in the distance matrix
  selected <- selection_matrix_from_clusters(clusters)
  N <- nrow(distances)
  for (i in 1:N) {
    exchange_partners <- get_exchange_partners(clusters, i, categories)
    ## Skip if there are zero exchange partners
    if (length(exchange_partners) == 0) {
      next
    }
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, length(exchange_partners))
    ## Swap item i with all legal exchange partners and store objectives
    for (j in seq_along(exchange_partners)) {
      comparison_objectives[j] <- update_objective_distance(
        distances = distances,
        selection_matrix = selected,
        i = i,
        j = exchange_partners[j],
        current_objective = best_total
      )
    }
    # Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # Which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      # Swap the elements in clustering vector
      clusters <- cluster_swap(clusters, i, swap)
      # Swap the elements in selection matrix
      selected <- swap_items(selected, i, swap)
      # Update best solution
      best_total <- best_this_round
    }
  }
  clusters
}

# Update objective value after swapping two items
#
# param distances: The N x N distance matrix
# param selection_matrix: An N x N matrix encoding the currently selected items
#                (generated via `selection_matrix_from_clusters`)
# param i, j: The items to be swapped
# param current_objective: The best found objective
# return: the objective value after swapping
update_objective_distance <- function(distances, selection_matrix, i, j, current_objective) {
  tmp_selection <- swap_items(selection_matrix, i, j)
  old_distances <- itemwise_distance_sum(distances, selection_matrix, i, j)
  new_distances <- itemwise_distance_sum(distances, tmp_selection, i, j)
  current_objective - old_distances + new_distances
}

# Create a boolean matrix for indexing in a distance matrix
# Inverse function for `clusters_from_selection_matrix`
#
# param clusters: A clustering vector with elements 1, ..., K indicating cluster membership
# return: A N x N matrix where TRUE in cell [i,j] indicates that elements
#     i and j are in the same (anti)cluster
selection_matrix_from_clusters <- function(clusters) {
  n <- length(clusters)
  K <- length(unique(clusters))
  selected <- matrix(FALSE, nrow = n, ncol = n)
  change_thing <- function(k) {
    in_same_cluster <- which(clusters == k)
    combis <- expand.grid(in_same_cluster, in_same_cluster)
    selected[as.matrix(combis)] <<- TRUE
  }
  lapply(1:K, change_thing)
  selected
}

# Convert a selection matrix to clustering vector
# Inverse function for `selection_matrix_from_clusters`
#
# param selection_matrix: A N x N matrix where TRUE in cell [i,j]
#     indicates that elements i and j are in the same (anti)cluster
# return: A clustering vector with elements 1, ..., K indicating cluster membership
clusters_from_selection_matrix <- function(selection_matrix) {
  N <- nrow(selection_matrix)
  clusters <- rep(NA, )
  # assign first cluster to set and all its connections
  clusters[1] <- 1
  clusters[which(selection_matrix[1, ] == TRUE)] <- 1
  # Now iterate over the remaining elements
  for (i in 2:N) {
    # test if item i is already in a cluster
    if (!is.na(clusters[i])) {
      next
    }
    next_cluster_index <- max(clusters, na.rm = TRUE) + 1
    clusters[i] <- next_cluster_index
    clusters[which(selection_matrix[i, ] == TRUE)] <- next_cluster_index
  }
  clusters
}

## Swap items in a boolean matrix for indexing.
swap_items <- function(selected, i, j) {
  tmp1 <- selected[i, ]
  tmp2 <- selected[, i]
  selected[i, ] <- selected[j, ]
  selected[, i] <- selected[, j]
  selected[j, ] <- tmp1
  selected[, j] <- tmp2
  selected[i, j] <- FALSE
  selected[j, i] <- FALSE
  selected
}

# (not that) important: distances has 0 in diagonal and upper triangular
itemwise_distance_sum <- function(distances, selected, i, j) {
  n <- nrow(distances)
  sum(
    distances[i, 1:(i-1)][selected[i, 1:(i-1)]],
    distances[i:n, i][selected[i:n, i]],
    distances[j, 1:(j-1)][selected[j, 1:(j-1)]],
    distances[j:n, j][selected[j:n, j]]
  )
}
