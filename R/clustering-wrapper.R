
#' Create K balanced clusters
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features.
#' @param K How many clusters should be created.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are assigned? Defaults to TRUE
#' @param method One of "heuristic" or "exact". When "exact" is chosen,
#'     the function minimizes the distance criterion (i.e., solves
#'     balanced k-cluster editing); when "heuristic"
#'     is chosen, the function tries to minimize the variance criterion
#'     (see \code{\link{anticlustering}} for details). Use
#'     \code{method = "exact"} only for rather small n.
#'
#' @return A vector representing the cluster affiliation of all elements.
#'     Each cluster has the the same size.
#'
#' @export
#'
#' @examples
#'
#' data(iris)
#' # Only use numeric attributes
#' clusters <- clustering(iris[, -5], K = 3)
#' # Compare feature means by anticluster
#' by(iris[, -5], clusters, function(x) round(colMeans(x), 2))
#' # Plot the anticlustering
#' par(mfrow = c(1, 2))
#' plot_clusters(iris[, 1:2], clusters)
#' plot_clusters(iris[, 3:4], clusters)
#' par(mfrow = c(1, 1))
#'
#' ## Exact balanced cluster editing method
#' # Create artifical data
#' n_features <- 2
#' n_elements <- 30
#' features <- matrix(rnorm(n_elements * n_features), ncol = n_features)
#' ac <- clustering(features, K = 2, method = "exact")
#' plot_clusters(features, ac, within_connection = TRUE)
#'

balanced_clustering <- function(features, K, standardize = TRUE,
                                method = "heuristic") {
  features <- as.matrix(features)
  if (standardize) {
    features <- scale(features)
  }

  if (method == "exact") {
    solver <- solver_available()
    if (solver == FALSE)
      stop("One of the packages 'Rglpk', 'gurobi', or 'Rcplex' ",
           "must be installed for an exact clustering")
    clusters <- equal_sized_cluster_editing(features, K, solver)
  } else if (method == "heuristic") {
    clusters <- equal_sized_kmeans(features, K)
  } else {
    stop("Argument `method` must be 'heuristic' or 'exact'.")
  }
  return(clusters)
}
