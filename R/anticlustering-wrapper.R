
#' Anticlustering
#'
#' Create equal sized groups of elements (anticlusters) that are as
#' similar as possible.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and/or lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param K How many anticlusters should be created.
#' @param objective The objective to be maximized. The option
#'     "distance" (default) is used to optimize the anticluster
#'     editing objective; the option "variance" is used to optimize
#'     the k-means anticlustering objective. See details.
#' @param method One of "heuristic" or "exact". See details.
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created. Defaults to \code{TRUE}. See
#'     details.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are created? Defaults to \code{TRUE}.
#'     Standardization is done using the function \code{\link{scale}}
#'     using the default settings (mean = 0, SD = 1). This argument only
#'     works in combination with the çode{features} argument, not with
#'     \code{distances}.
#' @param nrep The number of repetitions for the random sampling
#'     heuristic.  This argument only has an effect if \code{method}
#'     is \code{"heuristic"}. It does not have an effect if
#'     \code{method} is \code{"exact"}.
#'
#' @return A vector representing the anticluster affiliation.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#'
#' @export
#'
#' @details
#'
#' This function is used to solve balanced K anticlustering. That is, K
#' groups of equal size are created in such a way that similarity of all
#' groups is maximized. Set similarity is assessed using one of two
#' objective functions:
#'
#' - k-means *variance* objective, setting \code{objective = "variance"}
#'
#' - cluster editing *distance* objective, setting \code{objective = "distance"}
#'
#' The k-means objective maximizes the variance within anticlusters. The
#' cluster editing objective maximizes the sum of pairwise distances
#' within anticlusters. If the argument \code{features} is passed
#' together with \code{objective = "distance"} If another distance
#' measure is preferred, pass a self-computed dissimiliarity matrix via
#' the argument \code{distances}. The optimal cluster editing objective
#' can be found via integer linear programming; for the k-means
#' objective, there is only a heuristic option. Vary the parameter
#' \code{method} to select a "heuristic" or "exact" computation.
#' 
#' Both of these objectives are maximized to establish sets that are
#' similar; minimization of the same objectives creates a clustering,
#' i.e., establishes sets sets such that elements are as similar as
#' possible within a set and as different as possible between sets, see
#' \code{\link{balanced_clustering}}.
#'
#' To obtain an optimal solution for anticluster editing, a linear
#' programming solver must be installed and usable from R.  The
#' `anticlust` package supports the open source GNU linear programming
#' kit (called from the package \code{Rglpk}) and the commercial solvers
#' gurobi (called from the package \code{gurobi}) and IBM CPLEX (called
#' from the package \code{Rcplex}). A license is needed to use one of
#' the commercial solvers. The optimal solution is retrieved by setting
#' \code{objective = "distance"}, \code{preclustering = FALSE}, and
#' \code{method = "exact"}. Use this combination of arguments only for
#' small problem sizes (maybe <= 30 elements).
#'
#' To relax the optimality condition, it is possible to set
#' \code{preclustering = TRUE}. In this case, the anticluster editing
#' objective is still optimized using integer linear programming, but a
#' preprocessing forbids very similar elements to be assigned to the
#' same anticluster. This approach can be used to work on larger problem
#' instances and the solution is usually still optimal or very close to
#' optimal.
#'
#' If no exact solution is required or the problem size is too large for
#' integer linear programming, a heuristic method, based on repeated
#' random sampling, is available. Across a specified number of runs,
#' anticlusters are assigned randomly and the best assignment is
#' returned. This method works for both anticluster editing and k-means
#' anticlustering. The sampling approach may also incorporate a
#' preclustering that prevents grouping very similar elements into the
#' same anticluster; use \code{preclustering = TRUE} to activate this
#' option, which is also the default. It is suggested that the
#' preclustering condition is activated for the random sampling approach
#' because it usually improves the quality of the solution.
#'
#' @examples
#'
#' ## Use anticlustering on the iris data set with the default settings:
#' # (a) Optimizes the distance objective
#' # (b) Heuristic method: random sampling
#' # (c) 10,000 sampling repetitions
#' # (d) Preclustering is activated
#' # (e) Anticlustering uses standardized features (each feature has
#' #     mean = 0 and SD = 1)
#'
#' data(iris)
#' # Only use numeric attributes
#' anticlusters <- anticlustering(iris[, -5], K = 3)
#' # Compare feature means by anticluster
#' by(iris[, -5], anticlusters, function(x) round(colMeans(x), 2))
#' # Plot the anticlustering
#' par(mfrow = c(1, 2))
#' plot_clusters(iris[, 1:2], anticlusters)
#' plot_clusters(iris[, 3:4], anticlusters)
#'
#' ## Exact anticlustering
#' # Create artifical data
#' n_features <- 2
#' n_elements <- 20
#' K <- 2
#' features <- matrix(rnorm(n_elements * n_features), ncol = n_features)
#' anticlustering(features, K = K, method = "exact",
#'                preclustering = FALSE, standardize = FALSE)
#'
#' # Enable preclustering
#' anticlustering(features, K = K, method = "exact",
#'                preclustering = TRUE, standardize = FALSE)
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
#'

anticlustering <- function(features = NULL, distances = NULL,
                           K, objective = "distance",
                           method = "heuristic", preclustering = TRUE,
                           standardize = TRUE, nrep = 10000) {

  input_handling_anticlustering(features, distances, K, objective,
                                method, preclustering,
                                standardize, nrep)

  ## Standardize feature values (for each feature, mean = 0, sd = 1)?
  if (argument_exists(features)) {
    features <- as.matrix(features)
    if (standardize) {
      features <- scale(features)
    }
    distances <- as.matrix(dist(features))
  } else {
    distances <- as.matrix(as.dist(distances))
  }

  ## Exact method using ILP
  if (method == "exact") {
    return(exact_anticlustering(distances, K, solver_available(),
                                preclustering))
  }

  ## Heuristic method - possibly using preclustering
  preclusters <- NULL
  if (preclustering == TRUE) {
    n_preclusters <- nrow(features) / K
    if (objective == "distance" && K == 2) {
      preclusters <- greedy_matching(distances)
    } else if (objective == "distance" && K > 2) {
      preclusters <- greedy_balanced_k_clustering(distances, K)
    }
    else if (objective == "variance") {
      preclusters <- equal_sized_kmeans(features, n_preclusters)
    }
  }
  heuristic_anticlustering(features, K, preclusters,
                           objective, nrep = nrep, distances)
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
input_handling_anticlustering <- function(features, distances,
                                          K, objective, method,
                                          preclustering, standardize,
                                          nrep) {

  ## Validate feature input
  if (argument_exists(features)) {
    validate_input(features, "features", c("data.frame", "matrix", "numeric"))
    features <- as.matrix(features)
    validate_input(features, "features", objmode = "numeric")

    validate_input(K, "K", "numeric", len = 1,
                   greater_than = 1, must_be_integer = TRUE)
    if (nrow(features) %% K != 0) {
      stop("K must be a divider of the number of elements.")
    }
  }

  validate_input(nrep, "nrep", "numeric", len = 1, greater_than = 0,
                 must_be_integer = TRUE)
  validate_input(method, "method", len = 1,
                 input_set = c("exact", "heuristic"))
  validate_input(objective, "objective", len = 1,
                 input_set = c("variance", "distance"))

  validate_input(preclustering, "preclustering", "logical", len = 1,
                 input_set = c(TRUE, FALSE))
  validate_input(standardize, "standardize", "logical", len = 1,
                 input_set = c(TRUE, FALSE))

  if (method == "exact") {
    solver <- solver_available()
    if (solver == FALSE) {
        stop("An exact solution was requested, but none of the linear ",
             "programming packages 'Rglpk', 'gurobi', or 'Rcplex' is ",
             "installed. Try out method = 'heuristic' or install ",
             "a linear programming solver. E.g., install the GNU "
             "linear programming kit. ",
             "Visit http://gnuwin32.sourceforge.net/packages/glpk.htm ",
             "if you are using windows; , "
             "use homebrew to install it on mac (brew install glpk); ",
             "use the following command to install it on Ubuntu: ",
             "sudo apt install libglpk-dev. Then, install the package Rglpk using ,"
             "install.packages(Rglpk)")
    }
  }

  if (objective == "variance" && method == "exact") {
    stop("There is no exact method to maximize the variance criterion. ",
         "Use method = 'distance' instead")
  }

  if (!argument_exists(features) && !argument_exists(distances)) {
    stop("One of the arguments 'features' or 'distances' must be given.")
  }

  if (argument_exists(features) && argument_exists(distances)) {
    stop("Only pass one of the arguments 'features' or 'distances'.")
  }

  if (argument_exists(distances) && objective == "variance") {
    stop("'distances' cannot be used if 'objective' is 'variance'.")
  }

  return(invisible(NULL))
}
