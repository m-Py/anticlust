
#' Anticlustering
#'
#' Create groups of elements (anticlusters) that are as similar as possible.
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to elements and columns correspond to features.
#' @param n_anticlusters How many anticlusters should be created.
#' @param standardize Boolean - should the features be
#'     standardized before anticlusters are created? Defaults to TRUE
#' @param objective The objective to be maximized, either
#'     "distance" (default) or "variance". See Details
#' @param method One of "heuristic" or "exact". Currently, an exact
#'     solution can only be obtained when the `objective` is "distance".
#'
#' @return A vector representing anticluster affiliation
#'
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @details TODO: Write details
#'
#' @examples
#'
#' # TODO: Write examples
#'
anticlustering <- function(features, n_anticlusters, standardize = TRUE,
                           objective = "distance", method = "heuristic") {
  features <- as.matrix(features)
  if (standardize) {
    features <- scale(features)
  }

  if (objective == "distance") {
    heuristic <- ifelse(method == "exact", 1, 3)
    solver <- solver_available()
    if (method == "exact" & solver == FALSE)
      stop("One of the packages 'Rglpk', 'gurobi', or 'Rcplex' must be installed for an exact solution")
    anticlusters <- distance_anticlustering(features, n_anticlusters, solver, FALSE, heuristic)
  } else if (objective == "variance") {
    if (method == "exact")
      stop("There is no exact method for maximizing the variance criterion")
    n_preclusters <- nrow(features) / n_anticlusters
    preclusters <- equal_sized_kmeans(features, n_preclusters)
    anticlusters <- heuristic_anticlustering(features, preclusters, objective, nrep = 1000)
  } else {
    stop("Argument objective must be 'distance' or 'variance'.")
  }
  return(anticlusters)
}

## Check if a solver package can be used
solver_available <- function() {
  solvers <- c("Rcplex", "gurobi", "Rglpk")
  pcks <- rownames(installed.packages())
  solvers_available <- solvers %in% pcks
  if (sum(solvers_available) == 0) # no solver available
    return(FALSE)
  return(solvers[solvers_available][1]) # pick only one solver
}
