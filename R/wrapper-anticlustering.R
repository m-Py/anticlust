
#' Anticlustering
#'
#' Create sets of elements (anticlusters) that are as similar as 
#' possible to each other, by maximizing the heterogeneity within anticlusters.
#'
#' @param x The data input. Can be one of two structures: (1) A data matrix
#'     where rows correspond to elements and columns correspond to
#'     features (a single numeric feature can be passed as a vector). (2)
#'     An N x N matrix dissimilarity matrix; can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param K How many anticlusters should be created. Alternatively:
#'     A vector of length N where entry describes the initial grouping of 
#'     an input element.
#' @param objective The objective to be maximized. The option "distance"
#'     (default) maximizes the cluster editing objective function; the
#'     option "variance" maximizes the k-means objective function. See
#'     details.
#' @param method One of "exchange" (default) or "ilp".See
#'     details.
#' @param preclustering Boolean. Should a preclustering be conducted
#'     before anticlusters are created? Defaults to \code{FALSE}. See
#'     details.
#' @param categories A vector, data.frame or matrix representing one or
#'     several categorical constraints. See details.
#'
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and K) to each input element.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @details
#'
#' This function is used to solve »K anticlustering«. That is,
#' K groups are created in such a way that all groups are
#' as similar as possible. In the standard case, groups of equal
#' size are returned. Adjust the \code{K} argument to create groups
#' of different size. 
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
#' \strong{Heuristic anticlustering}
#'
#' In the default case, a heuristic method is employed for anticlustering: 
#' The exchange method (\code{method = "exchange"}): Building on an
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
#' specifically implemented for large data sets.
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
#' \strong{Categorical constraints}
#'
#' The argument \code{categories} may induce categorical constraints.
#' The grouping variables indicated by \code{categories} will be
#' balanced out across anticlusters. Currently, this functionality is
#' only available in combination with the exchange method, but not with 
#' the exact ILP approach. Note that
#' it is currently \strong{not} possible to apply preclustering constraints
#' and categorical constraints at the same time.
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
#' ## Optimize the cluster editing (distance) criterion using the exchange method
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 3,
#'   objective = "distance",
#'   method = "exchange"
#' )
#' # Compare feature means by anticluster
#' by(iris[, -5], anticlusters, function(x) round(colMeans(x), 2))
#' # Compare standard deviations by anticluster
#' by(iris[, -5], anticlusters, function(x) round(apply(x, 2, sd), 2))
#'
#'
#' ## Incorporate categorical restrictions:
#' anticlusters <- anticlustering(
#'   iris[, -5],
#'   K = 2,
#'   categories = iris[, 5]
#' )
#' table(iris[, 5], anticlusters)
#' 
#' 
#' ## Use preclustering and variance (k-means) criterion on large data sets
#' N <- 1000
#' K = 2
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' ac <- anticlustering(
#'   lds, 
#'   K = K,
#'   objective = "variance",
#'   preclustering = TRUE
#' )
#' 
#' # The following is equivalent to setting `preclustering = TRUE`:
#' cl <- balanced_clustering(lds, K = N / K)
#' ac <- anticlustering(
#'   lds, 
#'   K = K,
#'   objective = "variance",
#'   categories = cl
#' )
#'
#' @references
#'
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59-96.
#'
#' Papenberg, M., & Klau, G. W. (2019, October 30). Using anticlustering
#' to partition a stimulus pool into equivalent parts.
#' https://doi.org/10.31234/osf.io/3razc
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213-218.
#'

anticlustering <- function(x, K, objective = "distance", method = "exchange",
                           preclustering = FALSE, categories = NULL) {

  ## Get data into required format
  input_handling_anticlustering(x, K, objective, method, preclustering, categories)

  ## Only deal with 1 data object (either features or distances)
  data <- process_input(x)

  ## Exact method using ILP
  if (method == "ilp") {
    return(exact_anticlustering(
      data,
      K,
      solver_available(),
      preclustering)
    )
  }

  # Some preprocessing: get objective function, preclusters and categories:
  obj_function <- get_objective_function(data, objective, K)
  # Preclustering and categorical constraints are both processed in the
  # variable `categories` after this step:
  categories <- get_categorical_constraints(data, K, preclustering, categories)

  ## Redirect to fast exchange method for k-means exchange
  if (class(objective) != "function" && objective == "variance" &&
      method == "exchange" && sum(is.na(K)) == 0) {
    return(fast_anticlustering(data, K, Inf, categories))
  }

  ## Redirect to fast exchange method for anticluster editing
  if (class(objective) != "function" && objective == "distance" &&
      method == "exchange" && sum(is.na(K)) == 0) {
    return(fast_exchange_dist(data, K, categories))
  }

  ## General heuristic optimization:
  exchange_method(data, K, obj_function, categories)
}

# Function that processes input and returns the data set that the
# optimization is conducted on as matrix (for exchange method)
# Returned matrix either represents distances or features.
process_input <- function(data) {
  if (!is_distance_matrix(data)) {
    data <- as.matrix(data)
    return(data)
  }
  as.matrix(as.dist(data))
}

# Ensure that a distance matrix is passed
convert_to_distances <- function(data) {
  if (!is_distance_matrix(data)) {
    distances <- as.matrix(dist(data))
  } else {
    distances <- data
  }
  distances
}

# Determine the objective function needed for the input
# The function returns a function. It is ensured that the function
# removes missing values before computing the objective. Missing values
# mean that there are NAs in the clustering vector. 
get_objective_function <- function(data, objective, K) {
  if (class(objective) == "function") {
    obj_function <- objective
  } else {
    ## Determine how to compute objective, three cases:
    # 1. Distance objective, distances were passed
    # 2. Distance objective, features were passed
    # 3. Variance objective, features were passed
    if (objective == "distance" && is_distance_matrix(data)) {
      obj_function <- distance_objective_
    } else if (objective == "distance" && !is_distance_matrix(data)) {
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

  obj
}



# Determines if preclustering constraints or categorical constraints
# are present. Returns either of them or NULL. The input validation
# ensures that at most one of the constraints is present when this
# function is called.
get_categorical_constraints <- function(data, K, preclustering, categories) {
  if (preclustering == TRUE) {
    return(get_preclusters(data, K))
  }
  if (argument_exists(categories)) {
    return(merge_into_one_variable(categories))
  }
  NULL
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
  categories <- factor(do.call(paste0, as.list(categories)))
  # sort as numeric to get consistent return value
  order_cluster_vector(to_numeric(categories))
}

## function that computes preclusters
get_preclusters <- function(data, K) {
  if (length(K) > 1) {
    K <- length(unique(K))
  }
  N <- nrow(data)
  nn_centroid_clustering(data, K)
}
