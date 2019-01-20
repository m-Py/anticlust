
#' Based on a preclustering, create a heuristic item assignment
#'
#' @param features A data.frame representing the features that are used
#'   in the assignment.
#' @param clustering A vector representing the clustering of elements
#'
#' @return A vector representing the new group assignment.
#'
#' @export
#'
heuristic_item_assignment <- function(features, clustering, nrep = 100) {

  ## sort by group, later sort back by item and return group
  dat <- data.frame(group = clustering, features, item = 1:nrow(features))
  dat <- dat[order(dat$group), ]
  ## start optimizing
  best_obj <- -Inf
  n_items <- nrow(dat)
  n_groups <- n_items / length(unique(dat$group))
  ## sequentially try out random assignments that place pregrouped
  ## items in different groups
  for (i in 1:nrep) {
    group <- unlist(replicate(n_items / n_groups, sample(1:n_groups), simplify = FALSE))
    cur_obj <- obj_value(group, features)
    if (cur_obj > best_obj) {
      best_assignment <- group
      best_obj <- cur_obj
    }
  }
  dat$group <- best_assignment
  dat <- dat[order(dat$item), ]
  return(dat$group)
}

#' Objective value for an assignment of anticlusters
#'
#' Computes the objective value for the summed distances criterion
#'
#' @param assignment A vector representing the group assignment
#' @param features A data.frame representing the features that were used
#'   in the assignment.
#'
#' @return Scalar: the objective value.
#'
#' @export
#'

obj_value <- function(assignment, features) {
  ## determine distances within each group
  distances <- by(features, assignment, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}


#' Heuristic clustering algorithm creating equal-sized clusters
#'
#' Uses kmeans algorithm to initiate cluster centers and then
#' sequentially chooses a cluster center and assigns assigns the closest
#' item to it. This ensures that each cluster is filled with the same
#' number of items.
#'
#' @param items A vector, matrix or data.frame of data points.
#'     Rows correspond to items and columns correspond to features.
#' @param nclusters The number of clusters to be created. Must be a
#'     factor of the number of items.
#'
#' @return A vector representing the clustering
#'
#' @export
#'
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange
#' @importFrom stats kmeans
#'
equal_sized_clustering <- function(items, nclusters) {
  ## initialize cluster centers using kmeans
  centers <- stats::kmeans(items, nclusters)$centers
  ## determine distances between all items and all cluster centers:
  clust_dist <- dist_from_centers(items, centers)
  ## How often must I iterate to assign items to clusters? nitems /
  ## ncluster -> This is actually the number of groups in item
  ## assignment.
  nruns <- nrow(items) / nclusters
  ## storage of cluster assignments:
  assignments <- matrix(nrow = nrow(items) / nclusters, ncol = nclusters)
  ## iterate over runs (= number of groups in the item assignment
  ## problem)
  ## TODO: check if the two loops can be made faster
  for (i in 1:nruns) {
    ## iterate over clusters in random order
    for (j in sample(nclusters)) {
      # which cluster is addressed (using the random sequence):
      column <- paste0("cluster", j)
      ## which row has item that is closest:
      min_index <- which.min(clust_dist[, column])
      item <- clust_dist[min_index, ]$item
      ## remove the item from data.frame clust_dist:
      clust_dist <- clust_dist[-min_index, ]
      ## add item to cluster j
      assignments[i, j] <- item
    }
  }
  ## return in long format that is consistent with return of ILP
  assignments <- reshape2::melt(assignments)
  assignments <- assignments[, -1] # first column useless
  colnames(assignments) <- c("group", "item")
  assignments <- dplyr::arrange(assignments, item)
  return(assignments$group)
}

#' Edit distances of very close neighbours (similar items)
#'
#' Based on a preclustering, distances between items of the same cluster
#' are set to a large negative value. By considering a preclustering of
#' items, very similar items will not be assigned to the same group when
#' the fixed distance object is used to create the ILP formulation of
#' the item assignment instance.
#'
#' @param distances A distance object or matrix of
#'     between-item-distances.
#' @return A vector representing the group assignment; objects that
#'   are part of the same group as indicated by this vector are assigned
#'   a new distance.
#' @param value The value that is assigned to the fixed distances.
#'   Defaults to -1,000,000 currently.
#'
#' @return A distance object containing the fixed distances. For items
#'     that are part of the same cluster (as specified in the
#'     `data.frame` `assignment`), the between-item-distances are set to
#'     -1000000. This will have to be replaced by a theoretically-sound
#'     value; -1000000 is just a hack that will work in the present
#'     applications.
#'
#' @importFrom stats as.dist dist
#' @export

edit_distances <- function(distances, assignment, value = -1000000) {
  n_groups <- length(unique(assignment))
  distances <- as.matrix(distances)
  for (i in 1:n_groups) {
    items <- which(assignment == i) ## which items are in the i'th cluster
    ## two for-loops to tap all pairwise distances that may require
    ## fixing (a lot of unnecessary iterations, probably)
    for (j in 1:(length(items) - 1)) {
      for (t in 2:length(items)) {
        distances[items[j], items[t]] <- value
        distances[items[t], items[j]] <- value
      }
    }
  }
  return(as.dist(distances))
}
