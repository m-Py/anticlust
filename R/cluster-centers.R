
# Compute cluster centers
#
# @param features A data matrix of element features
# @param clusters A numeric vector indicating cluster membership of
#     each element
#
# @return A matrix of cluster centers. Rows represent clusters and
#   columns represent features
#

cluster_centers <- function(features, clusters) {
  centers <- by(features, clusters, colMeans, na.rm = TRUE)
  do.call(rbind, as.list(centers)) # as.list for the case of only one feature
}

# Determine distances of n data points to m cluster centers
#
#
# @param features A vector, matrix or data.frame of data points. If a
#     matrix or data.frame is passed, rows correspond to items and
#     columns to features.
# @param centers A matrix of cluster centers. Each row corresponds to a
#     cluster and each column corresponds to a feature (this format is,
#     for example, returned by the function `stats::kmeans` through the
#     element `centers`).
# @param squared Boolean - compute the squared euclidean distance?
#
# @return A data matrix; columns represent clusters
#     and contain the distance to the respective cluster for each item.
#
# @details
# This code was published in Leisch (2006).
#
# @references
# Leisch (2006). A Toolbox for K-Centroids Cluster Analysis. Computational
# Statistics and Data Analysis, 51(2), 526–544.
#
dist_from_centers <- function(features, centers, squared) {
  z <- matrix(0, nrow = nrow(features), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    if (squared)
      z[,k] <- colSums((t(features) - centers[k,])^2, na.rm = TRUE)
    else
      z[,k] <- sqrt( colSums((t(features) - centers[k,])^2, na.rm = TRUE) )
  }
  z
}

## A standard (or squared) euclidean distance between two data points
## Currently, this function is only used in the test files.
euc_dist <- function(x1, x2, squared = FALSE) {
  if (squared)
    return(sum((x1 - x2)^2))
  sqrt(sum((x1 - x2)^2))
}
