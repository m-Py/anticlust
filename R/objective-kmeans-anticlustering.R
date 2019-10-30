
#' Objective value for the variance criterion
#'
#' Check the objective value for a given clustering.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param clusters A vector representing (anti)clusters (e.g., returned
#'     by \code{\link{anticlustering}} or
#'     \code{\link{balanced_clustering}})
#'
#' @return The total within-cluster variance
#'
#' @details
#'
#' The variance objective is given by the sum of the squared
#' errors between cluster centers and individual data points. It is the
#' objective function used in k-means clustering, see
#' \code{\link{kmeans}}.
#'
#'
#' @references
#'
#' Jain, A. K. (2010). Data clustering: 50 years beyond k-means.
#' Pattern Recognition Letters, 31, 651–666.
#'
#' Papenberg, M., & Klau, G. W. (2019, October 30). Using anticlustering
#' to partition a stimulus pool into equivalent parts.
#' https://doi.org/10.31234/osf.io/3razc
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213–218.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' data(iris)
#' ## Clustering
#' clusters <- balanced_clustering(
#'   iris[, -5],
#'   K = 3
#' )
#' # This is low:
#' variance_objective(
#'   iris[, -5],
#'   clusters
#' )
#' ## Anticlustering
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   objective = "variance",
#'   nrep = 100
#' )
#' # This is higher:
#' variance_objective(
#'   iris[, -5],
#'   anticlusters
#' )
#'

variance_objective <- function(features, clusters) {
  validate_input(features, "features", c("data.frame", "matrix", "numeric"))
  features <- as.matrix(features)
  validate_input(features, "features", objmode = "numeric")
  validate_input(clusters, "anticlusters", class_string = c("numeric", "factor"))
  variance_objective_(clusters, features)
}

# Internal function - no input handling
variance_objective_ <- function(clusters, data) {
  ## 1. Compute cluster centers
  centers <- cluster_centers(data, clusters)
  ## 2. For each item, compute distance to each cluster center
  distances <- dist_from_centers(data, centers, squared = TRUE)
  ## 3. Use two-column matrix to select relevant distances
  distances <- distances[cbind(1:nrow(distances), clusters)]
  sum(distances)
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
