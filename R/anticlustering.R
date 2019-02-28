
#' Algorithms for anticlustering
#'
#' Create groups of elements (anticlusters) that are as similar as
#' possible.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features.
#'     A vector represents a single feature.
#' @param n_anticlusters How many anticlusters should be created.
#' @param objective The objective to be maximized, either "distance"
#'     (default) or "variance". See Details.
#' @param method One of "annealing", "sampling", or "exact". See details.
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created. Defaults to `TRUE`.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are created? Defaults to `TRUE`.
#'     Standardization is done using the function \code{\link{scale}}.
#' @param nrep The number of repetitions used in the heuristic methods
#'     "sampling" or "annealing".
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
#' criterion, use \code{preclustering = FALSE}, \code{method = "exact"}, and
#' \code{objective = "distance"}. Use method exact only for small problem
#' sizes (< 30 elements). To use the exact approach, a linear
#' programming solver must be installed and usable from R. The open
#' source GNU linear programming kit (called from the package \code{Rglpk})
#' or one of the commercial solvers gurobi (called from the package
#' \code{gurobi}) or IBM CPLEX (called from the package \code{Rcplex}) can be
#' used. A license is needed for the commercial solvers and one of the
#' interface packages must be installed.
#'
#' Two heuristic approaches are available, one based on repeated random
#' sampling and another based on simulated annealing. Both approaches
#' may rely on a preclustering that prevents grouping very similar elements
#' into the same anticluster if we set \code{preclustering = TRUE}.
#'
#' @examples
#'
#' # Use Iris data set with the default settings
#' data(iris)
#' anticlusters <- anticlustering(iris[, -5], n_anticlusters = 3)
#' # Compare feature means by anticluster
#' by(iris[, -5], anticlusters, colMeans)
#' # Plot the anticlustering
#' par(mfrow = c(1, 2))
#' plot_clusters(iris[, 1:2], anticlusters) # see overlap
#' plot_clusters(iris[, 3:4], anticlusters)
#'
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

anticlustering <- function(features, n_anticlusters, objective = "distance",
                           method = "annealing", preclustering = TRUE,
                           standardize = TRUE, nrep = 10000) {

  input_handling_anticlustering(features, n_anticlusters, objective,
                                method, preclustering,
                                standardize, nrep)

  ## Standardize feature values (for each feature, mean = 0, sd = 1)?
  features <- as.matrix(features) # if only one feature was passed
  if (standardize) {
    features <- scale(features)
  }

  ## Use exact method to solve distance anticlustering (using ILP)
  if (method == "exact") {
    anticlusters <- exact_anticlustering(features, n_anticlusters,
                                         solver_available(),
                                         preclustering)
  } else {
    ## Use heuristic approach.
    ## TODO: Preclustering not necessary if preclustering = FALSE
    n_preclusters <- nrow(features) / n_anticlusters
    preclusters   <- equal_sized_kmeans(features, n_preclusters)
    anticlusters  <- heuristic_anticlustering(features, preclusters,
                                              objective, nrep = nrep,
                                              method = method,
                                              preclustering = preclustering)
  }
  names(anticlusters) <- NULL
  return(anticlusters)
}


## A function for validating the arguments passed to `anticlustering`. This
## ensures that:
## (a) All arguments have correct type
## (b) Method "exact" can only be used with objective = "distance"
## (c) A solver package is install if method = "exact"
## (d) A legal number of anticlusters was requested
input_handling_anticlustering <- function(features, n_anticlusters, objective,
                                          method, preclustering, standardize, nrep) {

  ## Validate features as input
  validate_input(features, "features", c("data.frame", "matrix", "numeric"))
  features <- as.matrix(features)
  validate_input(features, "features", objmode = "numeric")

  validate_input(n_anticlusters, "n_anticlusters", "numeric", len = 1, greater_than = 1)
  if (nrow(features) %% n_anticlusters != 0) {
    stop("The number of anticlusters must be a divider of the number of elements.")
  }

  validate_input(nrep, "nrep", "numeric", len = 1, greater_than = 0,
                 must_be_integer = TRUE)
  validate_input(method, "method", c("character", "factor"), len = 1,
                 input_set = c("exact", "sampling",  "annealing"))
  validate_input(objective, "objective", c("character", "factor"), len = 1,
                 input_set = c("variance", "distance"))

  validate_input(preclustering, "preclustering", "logical", len = 1,
                 input_set = c(TRUE, FALSE))
  validate_input(standardize, "standardize", "logical", len = 1,
                 input_set = c(TRUE, FALSE))

  if (method == "exact") {
    solver <- solver_available()
    if (solver == FALSE) {
      stop("One of the packages 'Rglpk', 'gurobi', or 'Rcplex' ",
           "must be installed for an exact solution")
    }
  }

  if (objective == "variance" && method == "exact") {
    stop("There is no exact method for maximizing the variance criterion")
  }
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
