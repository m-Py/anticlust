
## Try implementing Spaeth's (1987) algorithm

# Uses an "exchange method":

# "That method, applied for (2), improves some (random) initial
# partition by successively moving on trial each object from
# its cluster to all the other ones, and by shifting it, if there is
# any reduction at all, to that one where the first term on the right
# side of (4) decreases most, otherwise taking the next object and
# finally passing through all the objects until no further improvement
# occurs.

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
  legal_number_of_clusters(features, clusters) # check input
  centers <- by(features, clusters, colMeans)
  do.call(rbind, as.list(centers)) # as.list for the case of only one feature
}

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
dist_from_centers <- function(features, centers) {
  features <- matrix(features)
  distances <- matrix(NA, nrow = nrow(features), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    distances[, k] <- sqrt(colSums((t(features) - centers[k, ])^2))
  }
  distances
}
