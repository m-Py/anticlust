

#' Heuristic algorithm creating equal-sized clusters
#'
#' Uses kmeans algorithm to initiate cluster centers, then sequentially
#' chooses a cluster center and assigns the closest element to it,
#' ensuring that each cluster is filled with the same number of items.
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to items and columns correspond to features.
#' @param nclusters The number of clusters to be created.
#'
#' @return A vector representing the clustering
#'
#' @examples
#'
#' @importFrom stats kmeans
#'
#' @noRd

equal_sized_kmeans <- function(features, nclusters) {
  ## kmeans does not deal with missing values
  if (any(is.na(features)))
    stop("The features include missing values.")
  ## initialize cluster centers using kmeans
  centers <- stats::kmeans(features, nclusters)$centers
  ## determine distances between all items and all cluster centers:
  clust_dist <- dist_from_centers(features, centers, squared = TRUE)
  assignments <- heuristic_cluster_assignment(clust_dist)
  return(clusters_to_long(assignments))
}

## Called from within `equal_sized_kmeans`. Assigns elements to
## clusters, sequentially fills elements to the closest cluster center,
## fills the same number of elements into each cluster
heuristic_cluster_assignment <- function(clust_dist) {
  ## Variables that encode relevant dimensions
  nclusters <- ncol(clust_dist)
  N <- nrow(clust_dist)
  elements_per_cluster <- N / nclusters
  ## Include element ID as column of input because rows are later
  ## removed from matrix
  clust_dist <- cbind(clust_dist, 1:N)
  colnames(clust_dist) <- c(paste0("cluster", 1:nclusters), "item")
  ## Storage of cluster assignments (Columns = clusters; cells indicate the
  ## index of the element that is assigned to the respective cluster)
  assignments <- matrix(nrow = elements_per_cluster, ncol = nclusters)
  for (i in 1:elements_per_cluster) {
    ## iterate over clusters in random order
    for (j in sample(nclusters)) {
      ## 1. Choose element that is closest to center of cluster j:
      min_index <- which.min(clust_dist[, j])
      element <- clust_dist[min_index, "item"]
      ## 2. Add this element to cluster j:
      assignments[i, j] <- element
      ## 3. Remove this element to prevent further assignments:
      clust_dist <- clust_dist[-min_index, , drop = FALSE] # drop = FALSE !
    }
  }
  return(assignments)
}

## Called from within `equal_sized_kmeans`. Turns a matrix of cluster
## assignments into long format.
clusters_to_long <- function(assignments) {
  P <- ncol(assignments) # number of clusters
  N <- P * nrow(assignments)
  N_PER_CLUSTER <- N / P
  long_clusters <- cbind(c(assignments), rep(1:P, each = N_PER_CLUSTER))
  ## Make in correct order:
  long_clusters <- sort_by_col(long_clusters, 1)
  return(long_clusters[, 2])
}


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
      z[,k] <- colSums((t(features) - centers[k,])^2)
    else
      z[,k] <- sqrt( colSums((t(features) - centers[k,])^2) )
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
