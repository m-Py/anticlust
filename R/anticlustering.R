
#' Algorithms for anticlustering
#'
#' Create groups of elements (anticlusters) that are as similar as
#' possible.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param n_anticlusters How many anticlusters should be created.
#' @param objective The objective to be maximized, either "distance"
#'     (default) or "variance". See details.
#' @param method One of "sampling", or "exact". See details.
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created. Defaults to \code{TRUE}. See
#'     details.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are created? Defaults to \code{TRUE}.
#'     Standardization is done using the function \code{\link{scale}}
#'     using the default settings.
#' @param nrep The number of repetitions used in the heuristic methods
#'     "sampling" or "annealing". This argument does not have an effect
#'     if the argument \code{method} is \code{"exact"}.
#'
#' @return A vector representing the anticluster affiliation.
#'
#' @importFrom utils installed.packages
#' @importFrom stats as.dist dist optim
#' @importFrom Matrix sparseMatrix
#'
#' @export
#'
#' @details
#'
#' Späth (1986) and Valev (1998) proposed to maximize the variance
#' criterion used in k-means clustering to establish similar groups in
#' the anticlustering application. In the \code{anticlust} package,
#' optimizing the variance criterion is accomplished by setting
#' \code{objective = "variance"}. The \code{anticlust} package also
#' introduces another objective function to the anticlustering
#' application that has been developed in the problem domain of cluster
#' editing and is based on a measure of the pairwise distances of data
#' points (Grötschel & Wakabayashi, 1989). In weighted cluster editing,
#' the optimal objective is found when the sum of within-cluster
#' distances is minimized; for the anticlustering application, the
#' distance objective is maximized instead.
#'
#' The \code{anticlust} uses integer linear programming to find optimal
#' objective for the distance criterion (Grötschel & Wakabayashi,
#' 1989). To obtain an optimal solution, a linear programming solver
#' must be installed on the system and be usable from R. The
#' \code{anticlust} package supports the open source GNU linear
#' programming kit (called from the package \code{Rglpk}) and the
#' commercial solvers gurobi (called from the package \code{gurobi}) and
#' IBM CPLEX (called from the package \code{Rcplex}). A license is
#' needed to use one of the commercial solvers.
#'
#' An optimal solution can only be obtained for distance anticlustering,
#' i.e., when setting \code{objective = "distance"}. To obtain the
#' optimal objective for the distance criterion, use the following
#' arguments: \code{preclustering = FALSE}, \code{method = "exact"}, and
#' \code{objective = "distance"}. Use this combination of arguments only
#' for small problem sizes (< 30 elements). To relax the optimality
#' condition, it is possible to set \code{preclustering = TRUE}. In this
#' case, the distance objective is still optimized using integer linear
#' programming, but a preprocessing forbids very similar elements to be
#' assigned to the same anticluster. This approach can be used to work
#' on larger problem instances and the solution is usually optimal or
#' very close to optimal.
#'
#' If no exact solution is required or the problem size is too large for
#' the exact approach, a heuristic method is available, based on
#' repeated random sampling. The approach may rely on a preclustering
#' that prevents grouping very similar elements into the same
#' anticluster if \code{preclustering = TRUE}. It is strongly suggested
#' that the preclustering condition is activated when using the heuristic
#' because the solution is usually better with preclustering.
#'
#' @examples
#'
#' ## Use anticlustering on the iris data set with the default settings:
#' # (a) Optimizes the distance objective
#' # (b) Heuristic method: random sampling
#' # (c) 10,000 sampling repetitions
#' # (d) Preclustering is activated
#' # (e) Anticlustering uses standardized features (Mean = 0, SD = 1)
#'
#' data(iris)
#' # Only use numeric attributes:
#' iris <- iris[, -5]
#' anticlusters <- anticlustering(iris, n_anticlusters = 3)
#' # Compare feature means by anticluster
#' by(iris[, -5], anticlusters, function(x) round(colMeans(x), 2))
#' # Plot the anticlustering
#' par(mfrow = c(1, 2))
#' plot_clusters(iris[, 1:2], anticlusters)
#' plot_clusters(iris[, 3:4], anticlusters)
#'
#' ## The exact approach
#' library("Rglpk") # package that calls the GNU linear programming kit
#' # Create artifical data
#' n_features <- 2
#' n_anticlusters <- 2
#' n_elements <- 20
#' features <- matrix(rnorm(n_elements * n_features), ncol = n_features)
#' ac <- anticlustering(features, n_anticlusters, method = "exact",
#'                      preclustering = FALSE, standardize = FALSE)
#' # Determine objective value (larger is better)
#' get_objective(features, ac, "distance")
#'
#' # Enable preclustering
#' ac_preclust <- anticlustering(features, n_anticlusters, method = "exact",
#'                               preclustering = TRUE, standardize = FALSE)
#' get_objective(features, ac_preclust, "distance")
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
#' Valev, V. (1998). Set partition principles revisited. In Joint IAPR
#' international workshops on statistical techniques in pattern
#' recognition (SPR) and structural and syntactic pattern recognition
#' (SSPR) (pp.  875–881).
#'

anticlustering <- function(features, n_anticlusters, objective = "distance",
                           method = "sampling", preclustering = TRUE,
                           standardize = TRUE, nrep = 10000) {

  input_handling_anticlustering(features, n_anticlusters, objective,
                                method, preclustering,
                                standardize, nrep)

  ## Standardize feature values (for each feature, mean = 0, sd = 1)?
  features <- as.matrix(features)
  if (standardize) {
    features <- scale(features)
  }

  ## Exact method using ILP
  if (method == "exact") {
    return(exact_anticlustering(features, n_anticlusters,
                                solver_available(), preclustering))
  }

  ## Heuristic method - possibly using preclustering
  preclusters <- NULL
  if (preclustering == TRUE) {
    n_preclusters <- nrow(features) / n_anticlusters
    preclusters   <- equal_sized_kmeans(features, n_preclusters)
  }
  heuristic_anticlustering(features, n_anticlusters, preclusters,
                           objective, nrep = nrep)

}


#' Validating the arguments passed to `anticlustering`
#'
#' This function ensures that:
#' (a) All arguments have correct type
#' (b) Method "exact" can only be used with objective = "distance"
#' (c) A solver package has to be installed if method = "exact"
#' (d) A legal number of anticlusters was requested
#'
#' Takes the same parameters as \code{anticlustering}
#'
#' @return NULL
#'
#' @noRd
input_handling_anticlustering <- function(features, n_anticlusters, objective,
                                          method, preclustering, standardize, nrep) {

  ## Validate feature input
  validate_input(features, "features", c("data.frame", "matrix", "numeric"))
  features <- as.matrix(features)
  validate_input(features, "features", objmode = "numeric")

  validate_input(n_anticlusters, "n_anticlusters", "numeric", len = 1,
                 greater_than = 1, must_be_integer = TRUE)
  if (nrow(features) %% n_anticlusters != 0) {
    stop("The number of anticlusters must be a divider of the number of elements.")
  }

  validate_input(nrep, "nrep", "numeric", len = 1, greater_than = 0,
                 must_be_integer = TRUE)
  validate_input(method, "method", len = 1,
                 input_set = c("exact", "sampling",  "annealing"))
  validate_input(objective, "objective", len = 1,
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
    stop("There is no exact method to maximize the variance criterion. ",
         "Use method = 'distance' instead")
  }
  return(invisible(NULL))
}


#' Check if a solver package can be used
#'
#' Has no parameters, checks in the installed packages from the user
#'
#' @return \code{FALSE} If no solver is available; A string identifying
#'   the solver if at least one is available ("Rcplex", "gurobi", "Rglpk")
#' @noRd
solver_available <- function() {
  solvers <- c("Rcplex", "gurobi", "Rglpk")
  pcks <- rownames(installed.packages())
  solvers_available <- solvers %in% pcks
  if (sum(solvers_available) == 0) # no solver available
    return(FALSE)
  return(solvers[solvers_available][1]) # pick only one solver
}
