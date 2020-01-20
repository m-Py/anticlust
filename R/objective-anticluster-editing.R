
#' Objective value for the cluster editing distance objective
#'
#' Check the objective value for a given clustering.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. A N x N matrix representing the
#'     pairwise dissimilarities between N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param clusters A vector representing (anti)clusters (e.g.,
#'     returned by \code{\link{anticlustering}} or
#'     \code{\link{balanced_clustering}}).
#'
#' @return The cluster editing objective
#'
#' @details
#'
#' The cluster editing objective objective is given by the sum of the
#' pairwise distances between elements within the same (anti)clusters.
#' When the argument \code{features} is passed, the Euclidean distance
#' is used.
#'
#'
#' @examples
#'
#' data(iris)
#' distances <- dist(iris[1:60, -5])
#' ## Clustering
#' clusters <- balanced_clustering(distances = distances, K = 3)
#' # This is low:
#' distance_objective(distances = distances, clusters = clusters)
#' ## Anticlustering
#' anticlusters <- anticlustering(distances = distances, K = 3)
#' # This is higher:
#' distance_objective(distances = distances, clusters = anticlusters)
#'
#'
#' # Illustrates the cluster editing objective as the sum of distances
#' # within groups (needs an integer linear programming solver!)
#' n_elements <- 12
#' features <- matrix(runif(n_elements * 2), ncol = 2)
#' n_groups <- 2
#' clusters <- balanced_clustering(features, K = n_groups, method = "ilp")
#' anticlusters <- anticlustering(features, K = n_groups, method = "ilp")
#' par(mfrow = c(1, 2))
#' plot_clusters(features, clusters, within_connection = TRUE,
#'               main = "Minimum within-group distances")
#' plot_clusters(features, anticlusters, within_connection = TRUE,
#'               main = "Maximum within-group distances")
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @references
#'
#' Papenberg, M., & Klau, G. W. (2019, October 30). Using anticlustering
#' to partition a stimulus pool into equivalent parts.
#' https://doi.org/10.31234/osf.io/3razc
#'

distance_objective <- function(features = NULL, distances = NULL,
                               clusters) {
  if (!argument_exists(features) && !argument_exists(distances)) {
    stop("One of the arguments 'features' or 'distances' must be given.")
  }

  if (argument_exists(features) && argument_exists(distances)) {
    stop("Only pass one of the arguments 'features' or 'distances'.")
  }
  validate_input(clusters, "anticlusters", class_string = c("numeric", "factor"))
  if (argument_exists(features)) {
    validate_input(features, "features", c("data.frame", "matrix", "numeric"))
    features <- as.matrix(features)
    validate_input(features, "features", objmode = "numeric")
    return(obj_value_distance(clusters, features))
  }
  validate_input(distances, "distances", c("matrix", "dist"))

  ## Compute the objective; the above only validates the input
  distances <- as.matrix(distances)
  distance_objective_(clusters, distances)
}


#' Internal function for distance objective via input through distances
#'
#' @param anticlusters A vector of n anticlusters
#' @param data A n x n matrix of distances
#'
#' @details
#' The second argument is named `data` to have a consistent
#' interface with the other objective value computation functions.
#'
#' @noRd

distance_objective_ <- function(anticlusters, data) {
  sums_within <- sapply(
    unique(anticlusters), 
    function(i) sum(as.dist(subset_data_matrix(data, anticlusters == i)))
  )
  return(sum(sums_within))
}


#' Objective value for the distance criterion
#'
#' @param anticlusters A vector representing the anticluster affiliation
#' @param data A data.frame, matrix or vector representing the
#'     features that are used in the assignment.
#'
#' @return Scalar, the total sum of within-cluster distances (based
#'     on the Euclidean distance).
#'
#' @details
#' The second argument is named `data` to have a consistent
#' interface with the other objective value computation functions.
#'
#' @importFrom stats dist
#'
#' @noRd

obj_value_distance <- function(anticlusters, data) {
  ## determine distances within each group
  distances <- by(data, anticlusters, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}
