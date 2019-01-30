
#' Anticlustering
#'
#' Create groups of elements (anticlusters) that are as similar as
#' possible.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features.
#' @param n_anticlusters How many anticlusters should be created.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are assigned? Defaults to FALSE
#' @param objective The objective to be maximized, either "distance"
#'     (default) or "variance". See Details.
#' @param method One of "annealing", "sampling", or "exact". See details.
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created. Defaults to TRUE and it is
#'     strongly suggested to keep it that way for the heuristic methods
#'     "annealing" and "sampling".
#'     `preclustering` = FALSE is mainly implemented to test against the
#'     option `preclustering` = TRUE.
#'
#' @return A vector representing anticluster affiliation
#'
#' @importFrom utils installed.packages
#' @importFrom stats as.dist dist optim
#' @importFrom Matrix sparseMatrix
#'
#' @export
#'
#' @details
#'
#' TODO: Some explanation of objective.
#'
#' An exact solution can only be obtained when the `objective` is
#' "distance". To obtain the optimal objective for the distance
#' criterion, use `preclustering` = FALSE, `method` = "exact", and
#' `objective` = "distance". Use method exact only for small problem
#' sizes (< 30 elements). To use the exact approach, a linear
#' programming solver must be installed and usable from R. The open
#' source GNU linear programming kit (called from the package `Rglpk`)
#' or one of the commercial solvers gurobi (called from the package
#' `gurobi`) or IBM CPLEX (called from the package `Rcplex`) can be
#' used. A license is needed for the commercial solvers and one of the
#' interface packages must be installed.
#'
#' Two heuristic approaches are available, one based on repeated random
#' sampling and another based on simulated annealing. Both approaches
#' rely on a preclustering that prevents grouping very similar elements
#' into the same anticluster. Method = "sampling" will be somewhat
#' faster, but method = "annealing" will usually return a slightly
#' better objective.
#'
#' @examples
#'
#' # Compare heuristic approaches
#' n_elements <- 200
#' features <- matrix(round(rnorm(n_elements * 2, 70, 15)), ncol = 2)
#' criterion <- "variance"
#' n_anticlusters <- 4
#' classes1 <- anticlustering(features, n_anticlusters, objective = criterion,
#'                            method = "sampling")
#' classes2 <- anticlustering(features, n_anticlusters, objective = criterion,
#'                            method = "annealing")
#' get_objective(features, classes1, objective = criterion)
#' get_objective(features, classes2, objective = criterion)
#'
#' # Try out exact approach, only works with distance objective
#' criterion <- "distance"
#' n_elements <- 20
#' features <- matrix(round(rnorm(n_elements * 2, 70, 15)), ncol = 2)
#' n_anticlusters <- 2
#' classes1 <- anticlustering(features, n_anticlusters, objective = criterion,
#'                            method = "sampling")
#' classes2 <- anticlustering(features, n_anticlusters, objective = criterion,
#'                            method = "annealing")
#' # exact variant
#' classes3 <- anticlustering(features, n_anticlusters, objective = criterion,
#'                            method = "exact", preclustering = FALSE)
#' get_objective(features, classes1, objective = criterion)
#' get_objective(features, classes2, objective = criterion)
#' get_objective(features, classes3, objective = criterion)
#'
#' @references
#'
#' M. Grötschel and Y. Wakabayashi, “A cutting plane algorithm for a
#' clustering problem,” Mathematical Programming, vol. 45, nos. 1-3, pp.
#' 59–96, 1989.
#'
#' H. Späth, “Anticlustering: Maximizing the variance criterion,”
#' Control and Cybernetics, vol. 15, no. 2, pp. 213-218, 1986.
#'

anticlustering <- function(features, n_anticlusters, standardize = FALSE,
                           objective = "distance", method = "annealing",
                           preclustering = TRUE) {

  if (!method %in% c("exact", "sampling",  "annealing"))
    stop("Method must be one of 'exact', 'sampling', or 'annealing'")
  if (!objective %in% c("variance", "distance"))
    stop("Argument objective must be 'distance' or 'variance'.")

  ## standardize feature values (for each feature, mean = 0, sd = 1)?
  if (standardize) {
    features <- scale(features)
  }

  ## 1. Distance objective is maximized
  if (objective == "distance" & method == "exact") {
    solver <- solver_available()
    if (method == "exact" & solver == FALSE)
      stop("One of the packages 'Rglpk', 'gurobi', or 'Rcplex' must be installed for an exact solution")
    heuristic  <- ifelse(preclustering == TRUE, 1, 0)
    anticlusters <- distance_anticlustering(features, n_anticlusters,
                                            solver, heuristic = heuristic)
  ## 2. Variance objective is maximized
  } else if (objective == "variance" & method == "exact") {
    stop("There is no exact method for maximizing the variance criterion")
  ## Heuristic approach for variance criterion:
  } else {
    n_preclusters <- nrow(features) / n_anticlusters
    preclusters <- equal_sized_kmeans(features, n_preclusters)
    anticlusters <- heuristic_anticlustering(features, preclusters,
                                             objective, nrep = 10000,
                                             method = method,
                                             preclustering = preclustering)
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
