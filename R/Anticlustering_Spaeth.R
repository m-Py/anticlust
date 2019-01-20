
## Try implementing Spaeth's (1987) algorithm

# Uses an "exchange method":

# "That method, applied for (2), improves some (random) initial
# partition by successively moving on trial each object from
# its cluster to all the other ones, and by shifting it, if there is
# any reduction at all, to that one where the first term on the right
# side of (4) decreases most, otherwise taking the next object and
# finally passing through all the objects until no further improvement
# occurs.

#' Compute objective value for variance criterion
#'
#' @param features A data.frame, matrix or vector representing the
#'   features that are used in the assignment.
#' @param anticlusters A vector representing the anticluster affiliation
#'
#' @return Scalar: the objective value
#'
#' @export
#'
obj_value_variance <- function(features, anticlusters) {
  ## 1. Compute cluster centers
  centers <- cluster_centers(features, anticlusters)
  ## 2. For each item, compute distance to each cluster center
  distances <- dist_from_centers(features, centers)
  ## 3. Determine distances to cluster centers within each anticluster
  summed_distances <- by(distances, INDICES = anticlusters, sum) # DOES NOT WORK AS EXPECTED!
  ## 4. Objective value is the sum of all distances per group
  return(sum(summed_distances))
}

#' Compute cluster centers
#'
#' @param features A data matrix of element features
#' @param clusters A numeric vector indicating cluster membership
#'   of each element
#'
#' @return A matrix of cluster centers. Rows represent clusters and
#'   columns represent features

cluster_centers <- function(features, clusters) {
  features <- as.matrix(features) #  if features is a vector
  centers <- by(features, clusters, colMeans)
  do.call(rbind, as.list(centers)) # as.list for the case of only one feature
}

# Determine distances of n data points to m cluster centers
#
# The data basis for the clustering algorithm implemented in function
# `equal_sized_clustering`.
#
# @param features A vector, matrix or data.frame of data points. If a
#     matrix or data.frame is passed, rows correspond to items and
#     columns to features.
# @param centers A matrix of cluster centers. Each row corresponds to a
#     cluster and each column corresponds to a feature (this format is,
#     for example, returned by the function `stats::kmeans` through the
#     element `centers`).
#
# @return A data matrix; columns represent clusters
#     and contain the distance to the respective cluster for each item.
#
dist_from_centers <- function(features, centers) {
  features <- as.matrix(features) # if points is only a vector
  ## store all distances from centers
  storage <- matrix(ncol = nrow(centers), nrow = nrow(features))
  ## determine distances from all cluster centers:
  for (i in 1:ncol(storage)) {
    storage[, i] <- dist_one_center(features, centers[i, ])
  }
  return(storage)
}

## compute the distances of a vector (or matrix or data.frame) of points
## to a cluster center. points has the same form as in the function
## TODO: this should be doable without for-loop to be faster!
dist_one_center <- function(points, center) {
  distances <- vector(length = nrow(points))
  for (i in 1:nrow(points)) {
    distances[i] <- euc_dist(points[i, ], center)
  }
  return(distances)
}


## A standard euclidian distance between two data points
euc_dist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}
