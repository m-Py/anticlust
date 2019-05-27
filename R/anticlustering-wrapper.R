
#' Anticlustering
#'
#' Create equal sized groups of elements (anticlusters) that are as
#' similar as possible.
#'
#' @param features A numeric vector, matrix or data.frame of data points.
#'     Rows correspond to elements and columns correspond to features. A
#'     vector represents a single numeric feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. A N x N matrix representing the
#'     pairwise dissimilarities between N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param K How many anticlusters should be created.
#' @param objective The objective to be maximized. The option
#'     "distance" (default) is used to optimize the anticluster
#'     editing objective; the option "variance" is used to optimize
#'     the k-means anticlustering objective. See details.
#' @param method One of "heuristic" or "ilp". See details.
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created? Defaults to \code{FALSE} See
#'     details.
#' @param standardize Boolean - should the features be standardized
#'     before anticlusters are created? Defaults to \code{FALSE}.
#'     Standardization is done using the function \code{\link{scale}}
#'     using the default settings (mean = 0, SD = 1). This argument
#'     only works when the input data is given via \code{features},
#'     not via \code{distances}.
#' @param nrep The number of repetitions for the random sampling
#'     heuristic. This argument only has an effect if \code{method}
#'     is \code{"heuristic"}. It does not have an effect if
#'     \code{method} is \code{"ilp"}.
#' @param categories A vector, data.frame or matrix representing
#'     one or several categorical constraints. These grouping
#'     variables are balanced out across anticlusters. Currently
#'     this functionality is only available in combination with the
#'     heuristic random sampling method and not with the ILP approach. If
#'     categorical contraints are employed, the value of the argument
#'     \code{preclustering} will be ignored (that is, there will be no
#'     preclustering).
#' @param parallelize Boolean. Indicates whether multiple processors should
#'     be used for the random sampling method.
#' @param seed A value to fixate the random seed when using the random
#'     sampling method. When \code{parallelize} is \code{TRUE}, using
#'     this argument is the only way to ensure reproducibility.
#'
#' @return A vector representing the anticluster affiliation of each
#'     input element.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#'
#' @export
#'
#' @details
#'
#' This function is used to solve »balanced K anticlustering«. That is,
#' K equal sized groups are created in such a way that similarity of
#' all groups is maximized. Set similarity is assessed using one of
#' two objective functions:
#'
#' - k-means *variance* objective, setting \code{objective =
#'   "variance"}
#'
#' - cluster editing *distance* objective, setting \code{objective =
#'   "distance"}
#'
#' The k-means objective maximizes the variance within anticlusters
#' (see \code{\link{variance_objective}}). The cluster editing
#' objective maximizes the sum of pairwise distances within
#' anticlusters (see \code{\link{distance_objective}}). Maximizing
#' either of these objectives is used to create similar groups;
#' minimization of the same objectives leads to a clustering, i.e.,
#' elements are as similar as possible within a set and as different
#' as possible between sets. Clustering could be done with the function
#' \code{\link{balanced_clustering}}.
#'
#' If the argument \code{features} is passed together with
#' \code{objective = "distance"}, the Euclidean distance between
#' elements is computed as the basic unit of the anticluster editing
#' objective. If another measure of dissimilarity is preferred, you
#' may pass a self-generated dissimiliarity matrix via the argument
#' \code{distances}.
#'
#' The optimal anticluster editing objective can be found via integer
#' linear programming. To this end, set \code{method = "ilp"}. To
#' obtain an optimal solution, a linear
#' programming solver must be installed and usable from R. The
#' `anticlust` package supports the open source GNU linear programming
#' kit (called from the package \code{Rglpk}) and the commercial
#' solvers gurobi (called from the package \code{gurobi}) and IBM
#' CPLEX (called from the package \code{Rcplex}). A license is needed
#' to use one of the commercial solvers. The optimal solution is
#' retrieved by setting \code{objective = "distance"},
#' \code{method = "ilp"} and \code{preclustering = FALSE}. Use this
#' combination of arguments only for small problem sizes (maybe <= 30
#' elements).
#'
#' To relax the optimality condition, it is possible to set the
#' argument \code{preclustering = TRUE}. In this case, the anticluster
#' editing objective is still optimized using integer linear
#' programming, but a preprocessing forbids very similar elements to
#' be assigned to the same anticluster. This approach can be used to
#' work on larger problem instances and the solution is usually still
#' optimal or very close to optimal.
#'
#' If no exact solution is required or the problem size is too large
#' for integer linear programming, a heuristic method is available via
#' setting \code{method = "heuristic"}. Note that this is the only
#' method when optimizing the variance criterion. The
#' heuristic employs repeated random sampling: across a specified
#' number of runs, anticlusters are assigned randomly and the best
#' assignment is returned. The sampling approach may also incorporate a
#' preclustering that prevents grouping very similar elements into the
#' same anticluster. Preclustering (for the exact and heuristic approach)
#' is performed by a call to \code{\link{balanced_clustering}} where the
#' argument \code{K} is set to the number of elements divided by the
#' number of anticlusters that are to be created (actually, this is not
#' exactly what happens internally, but it is equivalent).
#'
#' As the heuristic method relies on random sampling, the output will
#' vary between function calls. To make a computation results
#' reproducible, you can use the argument \code{seed} -- a specific
#' value will produce the same results when the function is called with
#' the same input parameters. Note that the
#' same seed will produce different results when the parameter
#' \code{parallelize} is varied. For the parallel computation, the
#' random seed is set using the function
#' \code{\link[parallel]{clusterSetRNGStream}}; otherwise, the function
#' \code{\link{set.seed}} is called.
#'
#'
#' @examples
#'
#' ## Use anticlustering on the iris data set. Create sets of plants
#' ## that are as similar as possible with regard to all four features
#' ## of the iris plants
#'
#' data(iris)
#' head(iris[, -5]) # these features are made similar across sets
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   nrep = 100 # increase for better results
#' )
#' # Compare feature means by anticluster
#' by(iris[, -5], anticlusters, function(x) round(colMeans(x), 2))
#'
#' # As the features are differently scaled, the solution might be improved
#' # by setting standardize = TRUE:
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   standardize = TRUE,
#'   nrep = 100
#' )
#'
#' # Optimize the variance criterion:
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   objective = "variance",
#'   nrep = 100
#' )
#'
#' # Increase nrep for better solutions; setting preclustering = TRUE
#' # sometimes also improves the solution
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   nrep = 100,
#'   preclustering = TRUE,
#'   objective = "variance"
#' )
#'
#' # Incorporate categorical restrictions:
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 2,
#'   categories = iris[, 5],
#'   nrep = 10
#' )
#' table(iris[, 5], anticlusters)
#'
#' @references
#'
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59-96.
#'
#' Böcker, S., Briesemeister, S., & Klau, G. W. (2011). Exact algorithms
#' for cluster editing: Evaluation and experiments. Algorithmica, 60,
#' 316-334.
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213-218.
#'
#'

anticlustering <- function(features = NULL, distances = NULL,
                           K, objective = "distance",
                           method = "heuristic", preclustering = FALSE,
                           standardize = FALSE, nrep = 10000,
                           categories = NULL, parallelize = FALSE,
                           seed = NULL) {

  input_handling_anticlustering(features, distances, K, objective,
                                method, preclustering,
                                standardize, nrep, categories,
                                parallelize, seed)

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
  if (method == "ilp") {
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
  if (argument_exists(categories)) {
    ## categorical constraints cannot be used with preclustering
    preclusters <- NULL
  }
  heuristic_anticlustering(features, K, preclusters,
                           objective, nrep = nrep, distances,
                           categories, parallelize, seed, ncores = NULL)
}


#' Validating the arguments passed to `anticlustering`
#'
#' This function ensures that:
#' (a) All arguments have correct type
#' (b) Method "ilp" can only be used with objective = "distance"
#' (c) A solver package has to be installed if method = "ilp"
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
                                          nrep, categories,
                                          parallelize, seed) {

  ## Validate feature input
  if (argument_exists(features)) {
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
                 input_set = c("ilp", "heuristic"))
  validate_input(objective, "objective", len = 1,
                 input_set = c("variance", "distance"))

  validate_input(preclustering, "preclustering", "logical", len = 1,
                 input_set = c(TRUE, FALSE))
  validate_input(standardize, "standardize", "logical", len = 1,
                 input_set = c(TRUE, FALSE))

  validate_input(parallelize, "parallelize", "logical", len = 1,
                 input_set = c(TRUE, FALSE))
  if (argument_exists(seed)) {
    validate_input(seed, "seed", "numeric", len = 1, not_na = TRUE)
  }

  if (method == "ilp") {
    solver <- solver_available()
    if (solver == FALSE) {
        stop("\n\nAn exact solution was requested, but none of the linear ",
             "programming \npackages 'Rglpk', 'gurobi', or 'Rcplex' is ",
             "available. \n\nTry `method = 'heuristic'` or install ",
             "a linear programming solver \nto obtain an exact solution. ",
             "For example, install the GNU linear \nprogramming kit: \n\n",
             "- On windows, visit ",
             "http://gnuwin32.sourceforge.net/packages/glpk.htm \n\n",
             "- Use homebrew to install it on mac, 'brew install glpk' \n\n",
             "- 'sudo apt install libglpk-dev' on Ubuntu ",
             "\n\nThen, install the Rglpk package via ",
             "`install.packages(Rglpk)`. \n\nOtherwise, you may obtain ",
             "a license for one of ",
             "the commercial solvers \ngurobi or IBM CPLEX (they are free ",
             "for academic use).")
    }
  }

  if (objective == "variance" && method == "ilp") {
    stop("You cannot use integer linear programming method to maximize the variance criterion. ",
         "Use objective = 'distance' or method = 'heuristic' instead")
  }

  if (!argument_exists(features) && !argument_exists(distances)) {
    stop("One of the arguments 'features' or 'distances' must be given.")
  }

  if (argument_exists(features) && argument_exists(distances)) {
    stop("Only pass one of the arguments 'features' or 'distances'.")
  }

  if (argument_exists(distances) && objective == "variance") {
    stop("The argument 'distances' cannot be used if the argument 'objective' is 'variance'.")
  }

  if (argument_exists(categories) && method == "ilp") {
    stop("The ILP method cannot incorporate categorical restrictions")
  }

  return(invisible(NULL))
}
