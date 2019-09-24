
#' Create K balanced clusters
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
#' @param K How many clusters should be created.
#' @param objective The objective to be minimized. The option
#'     "distance" (default) is used to minimize the cluster editing
#'     objective; the option "variance" is used to minimize the
#'     k-means objective. See details.
#' @param method One of "heuristic" or "ilp". See details.
#'
#' @return A vector representing the cluster affiliation of all elements.
#'
#' @details
#'
#' This function partitions a set of elements into K equal sized
#' clusters. The objective is one of the following:
#'
#' - k-means *variance* objective, setting \code{objective = "variance"}
#'
#' - cluster editing *distance* objective, setting \code{objective = "distance"}
#'
#' The k-means objective minimizes the variance within clusters. The
#' cluster editing objective minimizes the sum of pairwise distances
#' within clusters. If the argument \code{features} is passed together
#' with \code{objective = "distance"}, the Euclidean distance is
#' computed by default. If another distance measure is preferred, pass
#' a self-computed dissimiliarity matrix via the argument
#' \code{distances}. The optimal cluster editing objective can be
#' found via integer linear programming; for the k-means objective,
#' there is only a heuristic option. Vary the parameter \code{method}
#' to select a "heuristic" or exact (via "ilp") computation.
#'
#' To obtain an optimal solution for balanced cluster editing, a
#' linear programming solver must be installed and usable from R. The
#' `anticlust` package supports the open source GNU linear programming
#' kit (called from the package \code{Rglpk}) and the commercial
#' solvers gurobi (called from the package \code{gurobi}) and IBM
#' CPLEX (called from the package \code{Rcplex}). A license is needed
#' to use one of the commercial solvers. The optimal solution is
#' retrieved by setting \code{objective = "distance"} and \code{method
#' = "ilp"}.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' m.eik michalke \email{meik.michalke@@hhu.de}
#'
#' @examples
#'
#' data(iris)
#' # Only use numeric attributes
#' clusters <- balanced_clustering(iris[, -5], K = 3)
#' # Compare feature means by anticluster
#' by(iris[, -5], clusters, function(x) round(colMeans(x), 2))
#' # Plot the anticlustering
#' par(mfrow = c(1, 2))
#' plot_clusters(iris[, 1:2], clusters)
#' plot_clusters(iris[, 3:4], clusters)
#'
#'
#' ## Exact balanced cluster editing method
#' # Create artifical data
#' n_features <- 2
#' N <- 30
#' K <- 3
#' features <- matrix(rnorm(N * n_features), ncol = n_features)
#' ac1 <- balanced_clustering(features, K = K, method = "ilp")
#' ac2 <- balanced_clustering(features, K = K, method = "heuristic")
#'
#' ## Compare exact and heuristic balanced cluster editing
#' par(mfrow = c(1, 2))
#' plot_clusters(features, ac1, within_connection = TRUE,
#'               main = "optimal cluster editing", xlab = "", ylab = "")
#' plot_clusters(features, ac2, within_connection = TRUE,
#'               main = "heuristic cluster editing", xlab = "", ylab = "")
#'
#'
#' @references
#'
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59–96.
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213–218.
#'

balanced_clustering <- function(features = NULL, distances = NULL,
                                K, method = "heuristic") {

  input_handling_anticlustering(features, distances, K,
                                "distance", method, TRUE, 1,
                                NULL, NULL)

  if (argument_exists(features)) {
    features <- as.matrix(features)
    distances <- as.matrix(dist(features))
  } else if (argument_exists(distances)) {
    distances <- as.matrix(as.dist(distances))
  }

  if (method == "ilp") {
    return(balanced_cluster_editing(distances, K, solver_available()))
  }
  centroid_preclustering(features, distances, K = K)
}
