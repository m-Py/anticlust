
#' Cluster editing distance objective
#'
#' Compute the cluster editing objective for a given clustering.
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
#'     returned by \code{\link{anticlustering}}).
#'
#' @return The cluster editing objective
#'
#' @details
#'
#' The cluster editing objective objective is given by the sum of the
#' pairwise distances between elements within the same (anti)clusters.
#' When the input \code{x} is a feature matrix, the Euclidean distance
#' is computed as the basic distance unit.
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
  validate_input(clusters, "clusters", len = nrow(x), not_na = TRUE)
  distance_objective_(clusters, x)
}

# other order of the arguments, needed for some internal handling
distance_objective_ <- function(clusters, x) {
  sum(distance_objective_by_group(clusters, x))
}

# Compute distance objective by cluster
# param data: distance matrix or feature matrix
# param cl: cluster assignment
distance_objective_by_group <- function(cl, data) {
  if (is_distance_matrix(data)) {
    objectives <- sapply(
      sort(unique(cl)), 
      function(x) sum(as.dist(data[cl == x, cl == x]))
    )
  } else {
    objectives <- sapply(
      sort(unique(cl)), 
      function(x) sum(dist(data[cl == x, ]))
    )
  }
  objectives
}
