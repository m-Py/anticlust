
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
#' @param method One of "heuristic" or "ilp". See details.
#'
#' @return A vector representing the cluster affiliation of all elements.
#'
#' @details
#'
#' This function partitions a set of elements into K equal sized
#' clusters. The function offers two methods, a heuristic method and an
#' exact method. The heuristic (\code{method = "heuristic"}) computes
#' the centroid of all available elements and identifies the element
#' farthest to it (if the input is a dissimilarity matrix, the most
#' central element acts as the centroid). The farthest element is
#' clustered with its (N/K)-1 nearest neighbours. From the remaining
#' elements, the element farthest to the centroid is selected and again
#' clustered with its (N/K)-1 neighbours; the procedure is repeated
#' until all elements are part of a cluster.
#'
#' An exact method (\code{method = "ilp"}) can be used to solve cluster
#' editing optimally. The cluster editing objective minimizes the sum
#' of pairwise distances within clusters. If the argument
#' \code{features} is passed, the Euclidean distance is computed by
#' default as the basic unit of the cluster editing objective. If
#' another distance measure is preferred, users may pass a self-computed
#' dissimiliarity matrix via the argument \code{distances}. The optimal
#' cluster editing objective can be found via integer linear
#' programming. To obtain an optimal solution for balanced cluster
#' editing, a linear programming solver must be installed and usable
#' from R. The `anticlust` package supports the open source GNU linear
#' programming kit (called from the package \code{Rglpk}) and the
#' commercial solvers gurobi (called from the package \code{gurobi}) and
#' IBM CPLEX (called from the package \code{Rcplex}). A license is
#' needed to use one of the commercial solvers.
#'
#'
#' @source
#'
#' The heuristic method was developed and contributed by m.eik michalke.
#' The integer linear programming method was implemented by Martin
#' Papenberg.
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
  centroid_clustering(features, distances, K = K)
}
