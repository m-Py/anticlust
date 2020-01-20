
#' Objective value for the cluster editing distance objective
#'
#' Check the objective value for a given clustering.
#'
#' @param x The data input. Can be one of two structures: (1) A data matrix
#'     where rows correspond to elements and columns correspond to
#'     features (a single numeric feature can be passed as a vector). (2)
#'     An N x N matrix dissimilarity matrix; can be an object of class
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
#' clusters <- balanced_clustering(distances, K = 3)
#' # This is low:
#' distance_objective(distances, clusters)
#' ## Anticlustering
#' anticlusters <- anticlustering(distances, K = 3)
#' # This is higher:
#' distance_objective(distances, anticlusters)
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

distance_objective <- function(x, clusters) {
  x <- as.matrix(x)
  validate_input(x, "x", objmode = "numeric")

  if (!is_distance_matrix(x)) {
    return(obj_value_distance(clusters, x))
  }
  distances <- as.matrix(x)
  distance_objective_(clusters, x)
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
