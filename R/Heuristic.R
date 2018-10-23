
#' Based on a preclustering, this function creates a heuristic item assignment
#'
#' @param assignment A data.frame that represents the preclustering assingment.
#'   Usually, the return value of `ilp_to_groups` or `equal_sized_clustering`
#'   will be the input.
#'
#' @return The data.frame assigment that was passed, but now the column
#'   `group` contains the groups as determined for the item assignment
#'   problem.
#'
#'   @export
#'
heuristic_item_assignment <- function(assignment, nrep = 100) {
  assignment <- dplyr::arrange(assignment, assignment$group)
  best_obj <- -Inf
  n_items <- nrow(assignment)
  n_groups <- n_items / length(unique(assignment$group))

  ## sequentially try out random assignments that place pregrouped
  ## items in different groups. This procedure can be refined I guess
  for (i in 1:nrep) {
    assignment$group <- unlist(replicate(n_items / n_groups,
                                         sample(1:n_groups),
                                         simplify = FALSE))
    cur_obj <- obj_value(assignment)
    if (cur_obj > best_obj) {
      best_assignment <- assignment$group
      best_obj <- cur_obj
    }
  }
  assignment$group <- best_assignment
  assignment <- dplyr::arrange(assignment, assignment$item)
  return(assignment)
}

#' Get the objective function from an assignment object
#'
#' @param assignment A data.frame that represents a clustering
#'   assingment. Usually, the return value of `ilp_to_groups` or
#'   `equal_sized_clustering` will be the input.
#'
#' @return The objective value
#'
#' @export
#'

obj_value <- function(assignment) {
  ## determine distances within each group
  distances <- by(assignment[,3:ncol(assignment)], assignment$group, dist)
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
#' @param items A vector, matrix or data.frame of data points. If a
#'     matrix or data.frame is passed, rows correspond to items and
#'     columns correspond to features.
#' @param nclusters The number of clusters to be created. Must be a
#'     factor of the number of items.
#'
#' @return A data.frame containing one column of item ids (here, the id
#'     corresponds to the order of the items) and one column contains the
#'     group assignments of the items. The original items are also
#'     returned as columns of the data.frame.
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
  ## return feature values as well:
  assignments <- data.frame(assignments, items)
  return(assignments)
}


## A standard euclidian distance between two data points
euc_dist <- function(x1, x2)
  sqrt(sum((x1 - x2)^2))

# Determine distances of n data points to m cluster centers
#
# The data basis for the clustering algorithm implemented in function
# `equal_sized_clustering`.
#
# @param points A vector, matrix or data.frame of data points. If a
#     matrix or data.frame is passed, rows correspond to items and
#     columns to features.
# @param centers A matrix of cluster centers. Each row corresponds to a
#     cluster and each column corresponds to a feature (this format is,
#     for example, returned by the function `stats::kmeans` through the
#     element `centers`).
#
# @return A `data.frame`. The first column is named `item` and it just
#     contains a numbering of the rows (this is needed for the function
#     (`equal_sized_clustering`). The other columns represent clusters
#     and contain the distance to the respective column for each item.
#
dist_from_centers <- function(points, centers) {
  points <- data.frame(points) # if points is only a vector
  ## store all distances from centers
  storage <- matrix(ncol = nrow(centers), nrow = nrow(points))
  ## determine distances from all cluster centers:
  for (i in 1:ncol(storage)) {
    storage[, i] <- dist_one_center(points, centers[i, ])
  }
  ## add item number
  storage <- cbind(1:nrow(points), storage)
  ## add column names
  colnames(storage) <- c("item", paste0("cluster", 1:nrow(centers)))
  return(data.frame(storage))
}

## compute the distances of a vector (or matrix or data.frame) of points
## to a cluster center. points has the same form as in the function
## above.
dist_one_center <- function(points, center) {
  distances <- vector(length = nrow(points))
  for (i in 1:nrow(points)) {
    distances[i] <- euc_dist(points[i, ], center)
  }
  return(distances)
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
#' @param assignment A data.frame containing two columns: `group`
#'     contains the cluster assignment and `item` is the item id. The
#'     item id must correctly correspond to the rows in argument
#'     distances! The data.frame is usually created using exact cluster
#'     editing (i.e., using the function `solve_ilp` and setting the
#'     parameter `objective` to "min"), or by the function
#'     `equal_sized_clustering`.
#' @param value The value that is assigned to each fixed distance.
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
  n_groups <- length(unique(assignment$group))
  distances <- as.matrix(distances)
  for (i in 1:n_groups) {
    ith_group <- subset(assignment, group == i)
    items <- ith_group$item ## which items are in the i'th cluster
    ## two for-loops to tap all pairwise distances that may require
    ## fixing (a lot of unnecessary iterations, probably)
    for (j in 1:(length(items) - 1)) {
      for (t in 2:length(items)) {
        distances[items[j], items[t]] <- value # this is only a hack
        distances[items[t], items[j]] <- value
      }
    }
  }
  return(as.dist(distances))
}
