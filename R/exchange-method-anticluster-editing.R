
#' Solve anticluster editing using a fast exchange method
#'
#' @param distances An N x N dissimilarity matrix.
#' @param K The number of cluster or an initial cluster assignment
#' @param categories A vector representing categorical constraints
#' @param preclusters A vector representing preclustering constraints
#'
#' @return The anticluster assignment
#'
#' @noRd
#'
#'

fast_exchange_dist <- function(distances, K, categories) {
  distances[upper.tri(distances)] <- 0
  diag(distances) <- 0
  clusters <- initialize_clusters(distances, K, distance_objective_, categories)
  ## Matrix for indexing the elements in the distance matrix
  selected <- selection_matrix_from_clusters(clusters)
  N <- nrow(distances)
  best_total <- distance_objective_(clusters, distances)
  for (i in 1:N) {
    # cluster of current item
    group_i <- clusters[i]
    exchange_partners <- get_exchange_partners(clusters, i, categories)
    ## Do not use this item if there are zero exchange partners
    if (length(exchange_partners) == 0) {
      next
    }
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, length(exchange_partners))
    for (j in seq_along(exchange_partners)) {
      ## Swap item i with all legal exchange partners and check out objective
      comparison_objectives[j] <- update_objective_distance(
        distances = distances,
        selection_matrix = selected,
        i = i,
        j = exchange_partners[j],
        current_objective = best_total
      )
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      # swap the elements
      clusters[i] <- clusters[swap]
      clusters[swap] <- group_i
      ## also adjust the boolean selection matrix
      selected <- swap_items(selected, i, swap)
      # update best solution
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

## Initialize a boolean matrix for indexing in a distance matrix
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
