
#' Anticlustering
#'
#' Create equal-sized groups of elements (anticlusters) that are as
#' similar as possible.
#'
#' @param features A numeric vector, matrix or data.frame of data
#'     points.  Rows correspond to elements and columns correspond to
#'     features. A vector represents a single numeric feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. An N x N matrix representing the
#'     pairwise dissimilarities between N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param K How many anticlusters should be created.
#' @param objective The objective to be maximized. The option "distance"
#'     (default) maximizes the cluster editing objective function; the
#'     option "variance" maximizes the k-means objective function. See
#'     details.
#' @param method One of "exchange" (default), "sampling", or "ilp".  See
#'     details.
#' @param preclustering Boolean. Should a preclustering be conducted
#'     before anticlusters are created? Defaults to \code{FALSE}. See
#'     details.
#' @param standardize Boolean. Should the features be standardized
#'     before anticlusters are created? Defaults to \code{FALSE}.
#'     Standardization is done using the function \code{\link{scale}}
#'     using the default settings (mean = 0, SD = 1). This argument only
#'     works when the input data is given via \code{features}, not via
#'     \code{distances}.
#' @param nrep The number of repetitions for the random sampling
#'     method. Defaults to 10,000. This argument only has an effect if
#'     \code{method} is \code{"sampling"}.
#' @param categories A vector, data.frame or matrix representing one or
#'     several categorical constraints. See details.
#' @param parallelize Boolean. Indicates whether multiple processors
#'     should be used for the random sampling method. Defaults to
#'     \code{FALSE}.
#' @param seed A value to fixate the random seed when using the random
#'     sampling method. When \code{parallelize} is \code{TRUE}, using
#'     this argument is the only way to ensure reproducibility.
#'
#' @return A vector of length \code{nrow(features)} (or
#'     \code{nrow(distances)}) that assigns a group (i.e, a number
#'     between 1 and K) to each input element.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#'
#' @export
#'
#' @details
#'
#' This function is used to solve »balanced K anticlustering«. That is,
#' K groups of equal size are created in such a way that all groups are
#' as similar as possible. Set similarity is assessed using one of two
#' objective functions:
#'
#' - k-means *variance* objective, setting \code{objective = "variance"}
#'
#' - cluster editing *distance* objective, setting \code{objective =
#'   "distance"}
#'
#' The k-means objective maximizes the variance within
#' anticlusters---that is, the sum of the squared distances between each
#' element and its cluster center (see
#' \code{\link{variance_objective}}). The cluster editing objective
#' maximizes the sum of pairwise distances within anticlusters (see
#' \code{\link{distance_objective}}). Maximizing either clustering
#' objective will create create similar groups; minimization of the same
#' objectives leads to a clustering, i.e., elements are as similar as
#' possible within a set and as different as possible between
#' sets. (Such a clustering is also possible with the function
#' \code{\link{balanced_clustering}}.)
#'
#' If the argument \code{features} is passed together with
#' \code{objective = "distance"}, the Euclidean distance is computed as
#' the basic unit of the anticluster editing objective. If a different
#' measure of dissimilarity is preferred, you may pass a self-generated
#' dissimiliarity matrix via the argument \code{distances}.
#'
#' \strong{Exact anticlustering}
#'
#' The optimal anticluster editing objective can be found via integer
#' linear programming. To this end, set \code{method = "ilp"}. To obtain
#' an optimal solution, a linear programming solver must be installed
#' and usable from R. The \code{anticlust} package supports the open
#' source GNU linear programming kit (called from the package
#' \code{Rglpk}) and the commercial solvers gurobi (called from the
#' package \code{gurobi}) and IBM CPLEX (called from the package
#' \code{Rcplex}). A license is needed to use one of the commercial
#' solvers. The optimal solution is retrieved by setting \code{objective
#' = "distance"}, \code{method = "ilp"} and \code{preclustering =
#' FALSE}. Use this combination of arguments only for small problem
#' sizes (maybe <= 30 elements).
#'
#' To relax the optimality condition, it is possible to set the argument
#' \code{preclustering = TRUE}. In this case, the anticluster editing
#' objective is still optimized using integer linear programming, but a
#' preprocessing forbids very similar elements to be assigned to the
#' same anticluster. Thus, before the anticlustering objective is
#' optimized, a cluster analysis identifies small groups of similar
#' elements (pairs if K = 2, triplets if K = 3, and so forth). The
#' preclustering reduces the size of the solution space, making the ILP
#' approach applicable for larger problem instances. With preclustering,
#' optimality is no longer guaranteed, but the solution is usually
#' optimal or very close to optimal.
#'
#' The variance criterion cannot be solved to optimality using integer
#' linear programming. However, it is possible to employ the function
#' \code{\link{generate_partitions}} to obtain optimal solutions for
#' small problem instances.
#'
#' \strong{Heuristic anticlustering}
#'
#' In addition to the exact approach---that is only feasible for small
#' N---the function employs two heuristic approaches. The first option
#' is repeated random sampling (\code{method = "sampling"}): Across a
#' specified number of runs, each element is assigned to an anticluster
#' at random and the objective value associated with this assignment is
#' computed. In the end, the best assignment---the assignment that
#' maximized the objective function---is returned.  The second heuristic
#' is the exchange method (\code{method = "exchange"}): Building on an
#' initial random assignment, elements are swapped between anticlusters
#' in such a way that each swap improves set similarity by the largest
#' amount that is possible in each situation (cf. Späth, 1986). The
#' swapping procedure is repeated for each element; because each
#' possible swap is investigated for each element, the total number of
#' exchanges grows quadratically with input size, rendering the exchange
#' method unsuitable for large N. Setting \code{preclustering = TRUE}
#' will limit the legal exchange partners to very similar elements,
#' resulting in improved run time while preserving a rather good
#' solution. This option is recommended for larger N. For very large N,
#' check out the function \code{\link{fast_anticlustering}} that was
#' specifically implemented for large data sets, or use the random
#' sampling method.
#'
#' The random sampling approach may also incorporate a preclustering
#' that prevents grouping very similar elements into the same
#' anticluster (setting \code{preclustering = TRUE}). For large N (N >
#' 100), the preclustering tends to improve the quality of the solution
#' because the random sampling is conducted on a restricted solution
#' space that consists of rather good solutions.  Preclustering (for the
#' exact and heuristic approaches) is performed by a call to
#' \code{\link{balanced_clustering}} where the argument \code{K} is set
#' to the number of elements divided by the number of anticlusters that
#' are to be created (actually, this is not exactly what happens
#' internally, but it is equivalent).
#'
#' For the random sampling method, the output will vary between function
#' calls. To make a computation results reproducible, you can use the
#' argument \code{seed} -- a specific value will produce the same
#' results when the function is called with the same input
#' parameters. Note that the same seed will produce different results
#' when the parameter \code{parallelize} is varied. For the parallel
#' computation, the random seed is set using the function
#' \code{\link[parallel]{clusterSetRNGStream}}; otherwise, the function
#' \code{\link{set.seed}} is called.
#'
#' \strong{Recommendations}
#'
#' The following recommendations are provided on using this function:
#'
#' \enumerate{
#'  \item If it is desired that the feature means are as similar as possible
#'      between sets, select \code{objective = "variance"}
#'  \item If the average similarity between elements in different sets
#'      should be maximized (i.e., making sets as a whole similar to
#'      each other), select \code{objective = "distance"}
#'  \item When the objective is \code{"distance"} and an exact approach
#'      is infeasible: select \code{method = "exchange"}; it is also
#'      possible to activate \code{preclustering = TRUE}, which improves
#'      run time and often does not even impair quality of the solution.
#'  \item Generally, the exchange method is preferred over random
#'      sampling. Use random sampling only for large data sets or if
#'      a fast solution is needed. However, even in that case, using
#'      the exchange method with preclustering activated or using the function
#'      \code{\link{fast_anticlustering}} is usually preferred.
#' }
#'
#' \strong{Categorical constraints}
#'
#' The argument \code{categories} may induce categorical constraints.
#' The grouping variables indicated by \code{categories} will be
#' balanced out across anticlusters. Currently, this functionality is
#' only available in combination with the random sampling and exchange
#' method, but not with the exact ILP approach. Note that for the random
#' sampling method, it is \strong{not} possible to apply preclustering
#' constraints and categorical constraints at the same time. Instead,
#' the argument \code{preclustering} will be ignored (that is, there
#' will be no preclustering). The exchange method will try to adhere to
#' both preclustering and categorical constraints if \code{preclustering
#' = TRUE} while giving priority to the categorical constraints.
#'
#'
#' @seealso
#'
#' \code{\link{fast_anticlustering}}
#'
#' \code{\link{variance_objective}}
#'
#' \code{\link{distance_objective}}
#'
#' \code{\link{balanced_clustering}}
#'
#' \code{\link{generate_partitions}}
#'
#' @examples
#'
#' ## Use anticlustering on the iris data set. Create sets of plants
#' ## that are as similar as possible with regard to all four features
#' ## of the iris plants
#'
#' data(iris)
#' head(iris[, -5]) # these features are made similar across sets
#'
#'
#' ## Optimize the variance criterion using the exchange method
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   objective = "variance",
#'   method = "exchange"
#' )
#' # Compare feature means by anticluster
#' by(iris[, -5], anticlusters, function(x) round(colMeans(x), 2))
#' # Compare standard deviations by anticluster
#' by(iris[, -5], anticlusters, function(x) round(apply(x, 2, sd), 2))
#'
#'
#' ## Optimize the cluster editing objective while activating preclustering
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   method = "exchange",
#'   objective = "distance",
#'   preclustering = TRUE
#' )
#' by(iris[, -5], anticlusters, function(x) round(colMeans(x), 2))
#' # Anticluster editing (method = "distance") tends to make the
#' # standard deviations more similar, while k-means (method = "variance")
#' # tends to make the means more similar):
#' by(iris[, -5], anticlusters, function(x) round(apply(x, 2, sd), 2))
#'
#'
#'
#' ## Incorporate categorical restrictions:
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 2,
#'   categories = iris[, 5],
#'   method = "sampling",
#'   nrep = 1
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

anticlustering <- function(features = NULL, distances = NULL,
                            K, objective = "distance",
                            method = "exchange", preclustering = FALSE,
                            standardize = FALSE, nrep = 10000,
                            categories = NULL, parallelize = FALSE,
                            seed = NULL) {
  anticlustering_(features, distances, K, objective,
                  method, preclustering, standardize, nrep,
                  categories, parallelize,
                  seed, k_neighbours = Inf)
}

anticlustering_ <- function(features = NULL, distances = NULL,
                            K, objective = "distance",
                            method = "sampling", preclustering = FALSE,
                            standardize = FALSE, nrep = 10000,
                            categories = NULL, parallelize = FALSE,
                            seed = NULL, k_neighbours) {

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
    if (class(objective) != "function" && objective == "distance") {
      distances <- as.matrix(dist(features))
    }
  } else {
    distances <- as.matrix(as.dist(distances))
  }

  ## Exact method using ILP
  if (method == "ilp") {
    return(exact_anticlustering(distances, K, solver_available(),
                                preclustering))
  }
  ## Heuristic methods
  # Determine if "objective" was a function
  if (class(objective) != "function") {
    obj_function <- get_objective_function(features, distances, objective)
  } else {
    obj_function <- objective
  }
  if (is.logical(preclustering) && preclustering == TRUE) {
    preclusters <- get_preclusters(features, distances, K)
  } else if (is.logical(preclustering) && preclustering == FALSE) {
    preclusters <- NULL
  } else {
    preclusters <- preclustering # `preclustering` was already a preclustering vector
  }

  ## direct exchange method for k-means criterion to fast exchange for
  ## fast computation
  if (class(objective) != "function" &&
      is.null(preclusters) &&
      objective == "variance" &&
      method == "exchange") {
    method <- "fast-exchange"
  }

  categories <- merge_into_one_variable(categories) # may be NULL
  if (method == "sampling" || method == "heuristic") {
    if (argument_exists(categories)) {
      preclusters <- NULL
    }
    return(random_sampling(features, K, preclusters,
                           obj_function, nrep = nrep, distances,
                           categories, parallelize, seed,
                           ncores = NULL))
  } else if (method == "exchange") {
    return(exchange_method(features, distances, K, obj_function, categories, preclusters))
  } else if (method == "fast-exchange") {
    neighbours <- get_neighbours(features, k_neighbours, categories)
    clusters <- random_sampling(features, K, NULL, obj_function,
                                1, distances, categories, FALSE,
                                NULL, NULL)
    fast_exchange_(features, clusters, categories, neighbours)
  }
}

## Extracted function that computes the preclusters.
get_preclusters <- function(features, distances, K) {
  if (argument_exists(features)) {
    distances <- dist(features)
  }
  if (K == 2) {
    preclusters <- greedy_matching(distances)
  } else if (K > 2) {
    preclusters <- greedy_balanced_k_clustering(distances, K)
  }
  preclusters
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
    if (sum(!complete.cases(features)) >= 1) {
      warning("There are NAs in your data, take care!")
    }
    features <- as.matrix(features)
    validate_input(features, "features", objmode = "numeric")
    # allow that K is an initial assignment of elements to clusters
    if (length(K) == 1) {
      validate_input(K, "K", "numeric", len = 1,
                     greater_than = 1, must_be_integer = TRUE)
    } else {
      validate_input(K, "K", "numeric", len = nrow(features))
      if (method != "exchange") {
        stop("an initial cluster assignment only works with method = 'exchange'")
      }
    }
    if (length(K) == 1 && nrow(features) %% K != 0) {
      if (method == "ilp") {
        stop("K must be a divider of the number of elements with the ILP method. (Try out method = 'exchange' or method = 'sampling'.)")
      }
      if (is.logical(preclustering) && preclustering == TRUE) {
        stop("K must be a divider of the number of elements with preclustering. (Try out preclustering = FALSE.)")
      }
    }
  }

  validate_input(nrep, "nrep", "numeric", len = 1, greater_than = 0,
                 must_be_integer = TRUE)
  validate_input(method, "method", len = 1,
                 input_set = c("ilp", "sampling", "exchange", "heuristic", "fast-exchange"))

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
           "available. \n\nTry `method = 'sampling'`, `method = 'exchange'` or install ",
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

  if (class(objective) != "function" && objective == "variance" && method == "ilp") {
    stop("You cannot use integer linear programming method to maximize the variance criterion. ",
         "Use objective = 'distance', method = 'sampling', or method = 'exchange' instead")
  }

  if (!argument_exists(features) && !argument_exists(distances)) {
    stop("One of the arguments 'features' or 'distances' must be given.")
  }

  if (argument_exists(features) && argument_exists(distances)) {
    stop("Only pass one of the arguments 'features' or 'distances'.")
  }

  if (class(objective) != "function" && argument_exists(distances) && objective == "variance") {
    stop("The argument 'distances' cannot be used if the argument 'objective' is 'variance'.")
  }

  if (argument_exists(categories) && method == "ilp") {
    stop("The ILP method cannot incorporate categorical restrictions")
  }

  return(invisible(NULL))
}
