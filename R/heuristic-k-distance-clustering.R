
#' Greedy balanced k clustering
#'
#' Based on a distance input, creates equal sized clusters, trying
#' to minimize the sum of distances within clusters.
#'
#' @param distances A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}.
#' @param K The number of clusters.
#'
#' @return The clusters.
#'
#' @details
#'
#' A hierarchical cluster alogorithms that sequentielly add the element
#' to an existing cluster yielding the minimum sum of distances.
#'
#' @noRd
#'

greedy_balanced_k_clustering <- function(distances, K) {

  distances <- as.matrix(distances)
  N <- nrow(distances)

  ## Columns = anticlusters, rows = elements
  clusters <- matrix(nrow = N/K, ncol = K)
  ## Initialize first cluster elements
  first_elements <- init_clusters(distances, K)
  clusters[1, ] <- first_elements

  ## Set distances between elements that are already in a cluster to Inf
  up_to_now <- as.matrix(expand.grid(first_elements, first_elements))
  distances[up_to_now] <- Inf

  ## Fill the rest
  for (i in 2:nrow(clusters)) {
    for (k in 1:K) {
      ## Select items that are already in k'th cluster
      items_in_cluster <- clusters[, k]
      items_in_cluster <- items_in_cluster[!is.na(items_in_cluster)]
      ## Which is the new item that is inserted:
      insert <- best_cluster_fit(distances, items_in_cluster)
      ## Append new item to k'th anticluster
      clusters[i, k] <- insert
      ## Remove item from pool:
      complete <- c(clusters[!is.na(clusters)])
      up_to_now <- as.matrix(expand.grid(insert, c(complete)))
      distances[up_to_now] <- Inf
      distances[up_to_now[, 2:1]] <- Inf
    }
  }

  ## Convert matrix of clusters into vector
  groups <- rep(NA, N)
  for (k in 1:K) {
    groups[clusters[, k]] <- k
  }

  ## Assert that the output has the correct structure:
  if (any(table(groups) != table(groups)[2])) {
    stop("Something went wrong: not the same anticluster sizes")
  }
  if (any(is.na(groups))) {
    stop("Something went wrong: NA in output")
  }

  return(groups)
}

#' Determine the minimum (or maximum) sum of distances between
#' elements in a cluster and elements that are not yet in a cluster
#'
#' @param distances A distance matrix. It is in preprocessed form so
#'     that several entries will be Inf or -Inf. Is not of class `dist`,
#'     but is a `matrix`.
#' @param items_in_cluster The elements that are in the cluster under
#'     processing.
#' @param minimization Boolean; `TRUE` if we are looking for the closest
#'     element to a cluster. `FALSE` if we are looking for the element
#'     that is furthest away.
#'
#' @details
#'
#' The objective function that is optimized is the sum of distances
#' between one element and all elements that are already in a cluster.
#'
#' @return The element that is closest to (or furthest apart from)
#'     a cluster.
#'
#' @noRd
#'

best_cluster_fit <- function(distances, items_in_cluster,
                             minimization = TRUE) {
  ## Initialize variables encoding whether this is a min or max problem
  best_diff   <- ifelse(minimization, Inf, -Inf)
  best_item   <- NULL
  comparator  <- ifelse(minimization, `<`, `>`)

  ## Start comparing the fit of each item to the cluster
  for (i in 1:ncol(distances)) {
    ## Do not use columns of items that are already in the cluster
    if (i %in% items_in_cluster) {
      next
    }
    ## Sum of distances between current element and elements in cluster
    current_diff <- sum(distances[items_in_cluster, i])
    if (comparator(current_diff, best_diff)) {
      best_diff <- current_diff
      best_item <- i
    }
  }
  return(best_item)
}

#' Select the first K items that constitute the first element per cluster
#'
#' @param distances A distance matrix. It is in preprocessed form so
#'     that several entries will be Inf or -Inf. Is not of class `dist`,
#'     but is a `matrix`.
#' @param K The number of clusters
#'
#' @return A vector of length K representing the elements from different
#'     clusters.
#'
#' @noRd
#'
init_clusters <- function(distances, K) {
  clusters <- rep(NA, K)
  ## Pick the first two elements per anticluster
  clusters[1:2] <- which(distances == max(distances), arr.ind = TRUE)[1, ]
  if (K == 2) {
    return(clusters)
  }
  for (k in 3:K) {
    clusters[k] <- best_cluster_fit(distances, clusters[!is.na(clusters)], FALSE)
  }
  return(clusters)
}
