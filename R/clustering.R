
#' Create clusters of equal size
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features.
#' @param n_clusters How many anticlusters should be created.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are assigned? Defaults to TRUE
#' @param method One of "heuristic" or "exact". When "exact" is chosen,
#'     the function minimizes the distance criterion; when "heuristic"
#'     is chosen, the function tries to minimize the variance criterion
#'     (see details). Use method exact only for rather small problem
#'     sizes (< xyz elements).
#'
#' @return A vector representing the cluster affiliation of all elements.
#'
#' @importFrom stats kmeans
#'
#' @export
#'
#' @details TODO: Write details
#'
#' @examples
#'
#' # TODO: Write examples
#'
clustering <- function(features, n_clusters, standardize = TRUE, method = "heuristic") {
  features <- as.matrix(features)
  if (standardize) {
    features <- scale(features)
  }

  if (method == "exact") {
    solver <- solver_available()
    if (solver == FALSE)
      stop("One of the packages 'Rglpk', 'gurobi', or 'Rcplex' must be installed for an exact clustering")
    clusters <- equal_sized_cluster_editing(features, n_clusters, solver)
  } else if (method == "heuristic") {
    clusters <- equal_sized_kmeans(features, n_clusters)
  } else {
    stop("Argument `method` must be 'heuristic' or 'exact'.")
  }
  return(clusters)
}
