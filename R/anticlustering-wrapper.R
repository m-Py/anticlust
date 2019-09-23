
#' Anticlustering
#'
#' Create sets of elements (anticlusters) that are as similar as possible.
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
#' @param K How many anticlusters should be created. Alternatively:
#'     A vector of length N. Each entry in this vector
#'     describes the initial grouping of an input element (as a number
#'     between 1 and the number of groups).
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
#'     method. This argument only has an effect if \code{method} is
#'     \code{"sampling"}.
#' @param categories A vector, data.frame or matrix representing one or
#'     several categorical constraints. See details.
#' @param iv A vector, data.frame or matrix representing features that act
#'     as "independent variables", i.e., variables whose values
#'     are made as different as possible between sets (the opposite of
#'     the \code{features} argument). Cannot be used in conjunction with
#'     the argument \code{distances}.
#'
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and K) to each input element.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#'
#' @export
#'
#' @details
#'
#' This function is used to solve »K anticlustering«. That is,
#' K groups are created in such a way that all groups are
#' as similar as possible. In the standard case, groups of equal
#' size are returned. Adjust the \code{K} argument to create groups
#' of different size (see \code{\link{initialize_K}} for an example).
#'
#' Set similarity is assessed using one of two objective functions:
#'
#' - k-means *variance* objective, setting \code{objective = "variance"}
#'
#' - cluster editing *distance* objective, setting \code{objective =
#'   "distance"}
#'
#' The k-means objective measures the variance within
#' anticlusters---that is, the sum of the squared distances between each
#' element and its cluster center (see
#' \code{\link{variance_objective}}). The cluster editing objective
#' measures the sum of pairwise distances within anticlusters (see
#' \code{\link{distance_objective}}). Maximizing either of these
#' objectives will lead to similar groups (and this is what is
#' actually done when using this function). Minimization of the same
#' objectives would lead to a clustering, i.e., sets where elements are
#' similar within a set and different between sets.
#' (Such a clustering is also possible with the function
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
#' N---the function employs two heuristic approaches. One option
#' is repeated random sampling (\code{method = "sampling"}): Across a
#' specified number of runs, each element is assigned to an anticluster
#' at random and the objective value associated with this assignment is
#' computed. In the end, the best assignment---the assignment that
#' maximized the objective function---is returned.  The other option
#' is the exchange method (\code{method = "exchange"}): Building on an
#' initial random assignment, elements are swapped between anticlusters
#' in such a way that each swap improves set similarity by the largest
#' amount that is possible in a situation (cf. Späth, 1986). The
#' swapping procedure is repeated for each element; because each
#' possible swap is investigated for each element, the total number of
#' exchanges grows quadratically with input size, rendering the exchange
#' method unsuitable for large N. Setting \code{preclustering = TRUE}
#' will limit the legal exchange partners to very similar elements,
#' resulting in improved run time while preserving a rather good
#' solution. This option is recommended for larger N. For very large N,
#' check out the function \code{\link{fast_anticlustering}} that was
#' specifically implemented for large data sets (or use the random
#' sampling method with few repetitions).
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
#' \strong{Min-max anticlustering}
#'
#' It is possible to conduct »min-max anticlustering«. That is: For
#' some variables, similarity is maximized between sets, whereas for
#' other variables, similarity is minimized. Use the argument \code{iv}
#' to define variables whose values should be dissimilar between sets;
#' \code{iv} is used the same way as \code{features} whose values are
#' made similar between sets.
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
#' @seealso
#'
#' \code{\link{fast_anticlustering}}
#'
#' \code{\link{initialize_K}}
#'
#' \code{\link{variance_objective}}
#'
#' \code{\link{distance_objective}}
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
                           standardize = FALSE, nrep = 10,
                           categories = NULL, iv = NULL) {

  input_handling_anticlustering(features, distances, K, objective,
                                method, preclustering, standardize,
                                nrep, categories, iv)

  ## Exact method using ILP
  if (method == "ilp") {
    return(exact_anticlustering(
      features,
      distances,
      K,
      solver_available(),
      preclustering)
    )
  }

  ## Get data into required format and get objective function:
  categories <- merge_into_one_variable(categories) # may be NULL
  data <- process_input(features, distances, standardize, objective, method)
  obj_function <- get_objective_function(features, distances, objective, K, iv)
  preclusters <- get_preclusters(features, distances, K, preclustering)

  ## Redirect to fast exchange method for k-means exchange (and no preclustering)
  if (class(objective) != "function" && is.null(preclusters) &&
      objective == "variance" &&  method == "exchange" &&
      sum(is.na(K)) == 0 && !argument_exists(iv)) {
    return(fast_anticlustering(features, K, Inf, categories))
  }

  ## Redirect to fast exchange method for anticluster editing
  if (class(objective) != "function" && objective == "distance"
      && method == "exchange"  && !argument_exists(iv)) {
    return(fast_exchange_dist(data, K, categories, preclusters))
  }

  ## Start heuristic optimization:
  heuristic_anticlustering(data, K, obj_function,
                           method, preclusters, nrep,
                           categories)
}

## Function that processes input and returns the data set that the
## optimization is conducted on (for exchange and sampling methods)
process_input <- function(features, distances, standardize, objective, method) {
  if (argument_exists(features)) {
    data <- as.matrix(features)
    if (standardize) {
      data <- scale(data)
    }
    ## Why is the data only coverted to distances for exchange method?
    ## -> random sampling uses `obj_value_distance` which uses features
    ## as input. Exchange method operates on distances.
    if (class(objective) != "function" && objective == "distance" && method == "exchange") {
      data <- as.matrix(dist(features))
    }
    return(data)
  }
  as.matrix(as.dist(distances))
}

#' Determine the objective function needed for the input
#'
#' @noRd
get_objective_function <- function(features, distances, objective, K, iv) {
  if (class(objective) == "function") {
    obj_function <- objective
  } else {
    ## What was the input: features or distances
    use_distances <- FALSE
    if (argument_exists(distances)) {
      use_distances <- TRUE
    }
    ## Determine how to compute objective, three cases:
    # 1. Distance objective, distances were passed
    # 2. Distance objective, features were passed
    # 3. Variance objective, features were passed
    if (objective == "distance" && use_distances == TRUE) {
      obj_function <- distance_objective_
    } else if (objective == "distance" && use_distances == FALSE) {
      obj_function <- obj_value_distance
    } else {
      obj_function <- variance_objective_
    }
  }

  ## Handle NA in initial cluster assignment
  if (sum(is.na(K)) > 0) {
    obj <- function(clusters, data) {
      data <- data[!is.na(clusters), , drop = FALSE]
      clusters <- clusters[!is.na(clusters)]
      obj_function(clusters, data)
    }
  } else {
    obj <- obj_function
  }

  ## Min-max application: independent variable is present
  if (argument_exists(iv)) { # make means in independent variable different
    iv <- as.matrix(iv)
    obj2 <- function(clusters, data) {
      obj(clusters, data) + obj(clusters, iv) * (-1)
    }
  } else {
    obj2 <- obj
  }
  obj2
}


#' Merge several grouping variable into one
#'
#' @param categories A vector, data.frame or matrix that represents
#'     one or several categorical constraints.
#'
#' @return A vector representing the group membership (or the combination
#'     of group memberships) as one variable
#'
#' @noRd
#'

merge_into_one_variable <- function(categories) {
  if (is.null(categories)) {
    return(NULL)
  }
  categories <- data.frame(categories)
  factor(do.call(paste0, as.list(categories)))
}


## function that computes preclusters
get_preclusters <- function(features, distances, K, preclustering) {
  if (length(K) > 1) {
    K <- length(unique(K))
  }
  ## Get precluster, three cases are possible
  # (a) Preclusters may be NULL (if preclustering == FALSE)
  # (b) may need to be computed (if preclustering == TRUE)
  # (c) it was already passed by the user (preclustering = a vector)
  if (is.logical(preclustering) && preclustering == FALSE) {
    return(NULL)
  }
  if (is.logical(preclustering) && preclustering == TRUE) {
    if (argument_exists(features)) {
      distances <- dist(features)
    }
    N <- nrow(as.matrix(distances))
    if (K == 2) {
      preclusters <- greedy_matching(distances)
    } else if (K > 2) {
      preclusters <- greedy_balanced_k_clustering(distances, N / K)
    }
  } else {
    preclusters <- preclustering
  }
  preclusters
}

# Direct to exchange method or sampling
heuristic_anticlustering <- function(data, K, obj_function,
                                     method, preclusters, nrep,
                                     categories) {
  if (method == "sampling" || method == "heuristic") {
    # For random sampling, we cannot apply preclustering and categorical
    # constraints at the same time
    if (argument_exists(categories)) {
      preclusters <- NULL
    }
    return(random_sampling(data, K, preclusters,
                           obj_function, nrep, categories))
  } else if (method == "exchange") {
    return(exchange_method(data, K, obj_function, categories, preclusters))
  }
}
