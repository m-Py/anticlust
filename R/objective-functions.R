
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
#' @param standardize Boolean - should the features be standardized
#'     before the objective is computed? Defaults to \code{FALSE}.
#'     Standardization is done using the function \code{\link{scale}}
#'     using the default settings (mean = 0, SD = 1).
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
#' @note
#'
#' When using this function to check the results of
#' \code{\link{anticlustering}} or \code{\link{balanced_clustering}},
#' make sure that the \code{standardization} argument has the same value
#' when creating (anti)clusters and when calling \code{variance_objective}.
#'
#' @references
#'
#' Jain, A. K. (2010). Data clustering: 50 years beyond k-means.
#' Pattern Recognition Letters, 31, 651–666.
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213–218.
#'
#' @export
#'
#' @examples
#'
#' data(iris)
#' ## Clustering
#' clusters <- balanced_clustering(iris[, -5], K = 3, standardize = TRUE)
#' # This is low:
#' variance_objective(iris[, -5], clusters, standardize = TRUE)
#' ## Anticlustering
#' anticlusters <- anticlustering(iris[, -5], K = 3, standardize = TRUE)
#' # This is higher:
#' variance_objective(iris[, -5], anticlusters, standardize = TRUE)
#'

variance_objective <- function(features, clusters,
                               standardize = FALSE) {

  validate_input(features, "features", c("data.frame", "matrix", "numeric"))
  features <- as.matrix(features)
  validate_input(features, "features", objmode = "numeric")
  validate_input(clusters, "anticlusters", class_string = c("numeric", "factor"))
  validate_input(standardize, "standardize", input_set = c(TRUE, FALSE),
                 not_na = TRUE, len = 1)

  if (standardize) {
    features <- scale(features)
  }

  ## 1. Compute cluster centers
  centers <- cluster_centers(features, clusters)
  ## 2. For each item, compute distance to each cluster center
  distances <- dist_from_centers(features, centers, squared = TRUE)
  ## 3. Use two-column matrix to select relevant distances
  distances <- distances[cbind(1:nrow(distances), clusters)]
  return(sum(distances))
}

#' Objective value for the cluster editing distance objective
#'
#' Check the objective value for a given clustering.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and/or lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param clusters A vector representing (anti)clusters (e.g.,
#'     returned by \code{\link{anticlustering}} or
#'     \code{\link{balanced_clustering}}).
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are created? Defaults to \code{FALSE}.
#'     Standardization is done using the function \code{\link{scale}}
#'     using the default settings (mean = 0, SD = 1). This argument
#'     only works in combination with the \code{features} argument,
#'     not with \code{distances}.
#'
#' @return The cluster editing objective
#'
#' @details
#'
#' The cluster editing objective objective is given by the sum of the
#' pairwise distances between elements within the same (anti)clusters.
#'
#' @note
#'
#' When using this function to check the results of
#' \code{\link{anticlustering}} or \code{\link{balanced_clustering}},
#' make sure that the \code{standardization} argument has the same value
#' when creating (anti)clusters and when calling \code{variance_objective}
#' (at least if \code{features} is used as input).
#'
#' @examples
#'
#' data(iris)
#' distances <- dist(iris[, -5])
#' ## Clustering
#' clusters <- balanced_clustering(distances = distances, K = 3)
#' # This is low:
#' distance_objective(distances = distances, clusters = clusters)
#' ## Anticlustering
#' anticlusters <- anticlustering(distances = distances, K = 3)
#' # This is higher:
#' distance_objective(distances = distances, clusters = anticlusters)
#'
#' ## We can also use the features as input, using a standardization
#' clusters <- balanced_clustering(iris[, -5], K = 3, standardize = TRUE)
#' distance_objective(iris[, -5], clusters = clusters, standardize = TRUE)
#'
#' anticlusters <- anticlustering(iris[, -5], K = 3, standardize = TRUE)
#' distance_objective(iris[, -5], clusters = anticlusters, standardize = TRUE)
#'
#' # Illustrates the cluster editing objective as the sum of distances
#' # within groups:
#' n_elements <- 12
#' features <- matrix(runif(n_elements * 2), ncol = 2)
#' n_groups <- 3
#' clusters <- balanced_clustering(features, K = n_groups, method = "exact")
#' anticlusters <- anticlustering(features, K = n_groups, method = "exact")
#' par(mfrow = c(1, 2))
#' plot_clusters(features, clusters, within_connection = TRUE)
#' plot_clusters(features, anticlusters, within_connection = TRUE)
#'
#' @export
#'
#' @references
#'
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59–96.
#'

distance_objective <- function(features = NULL, distances = NULL,
                               clusters, standardize = FALSE) {
  if (!argument_exists(features) && !argument_exists(distances)) {
    stop("One of the arguments 'features' or 'distances' must be given.")
  }

  if (argument_exists(features) && argument_exists(distances)) {
    stop("Only pass one of the arguments 'features' or 'distances'.")
  }
  validate_input(clusters, "anticlusters", class_string = c("numeric", "factor"))
  validate_input(standardize, "standardize", input_set = c(TRUE, FALSE),
                 not_na = TRUE, len = 1)
  if (argument_exists(features)) {
    validate_input(features, "features", c("data.frame", "matrix", "numeric"))
    features <- as.matrix(features)
    validate_input(features, "features", objmode = "numeric")
    if (standardize) {
      features <- scale(features)
    }
    return(obj_value_distance(features, clusters))
  }
  validate_input(distances, "distances", c("matrix", "dist"))
  distance_objective_(distances, clusters)
}


#' Internal function for distance objective via input through distances
#' @noRd
distance_objective_ <- function(distances, anticlusters) {
  K <- length(unique(anticlusters))
  distances <- as.matrix(distances)
  sums_within <- rep(NA, K)
  for (k in 1:K) {
    elements <- which(anticlusters == k)
    selection <- unique_combinations(elements)
    sums_within[k] <- sum(distances[selection])
  }
  return(sum(sums_within))
}

## This variant to create unique combinations is *much* faster
## (especially for larger N) than using utils::combn even though it
## seems that a lot of unnecessary combinations are generated using
## expand.grid
unique_combinations <- function(elements) {
  selection <- as.matrix(expand.grid(elements, elements))
  selection[selection[, 1] < selection[, 2], ]
}


#' Objective value for the distance criterion
#'
#' Currently I am using `distance_objective_`, even though this function
#' may sometimes be faster.
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used in the assignment.
#' @param anticlusters A vector representing the anticluster affiliation
#'
#' @return Scalar, the total sum of within-cluster distances (based
#'     on the Euclidean distance).
#'
#' @importFrom stats dist
#'
#' @noRd

obj_value_distance <- function(features, anticlusters) {
  ## determine distances within each group
  distances <- by(features, anticlusters, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}
