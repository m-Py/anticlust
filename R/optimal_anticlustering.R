#' Optimal ("exact") algorithms for anticlustering
#' 
#' Wrapper function that gives access to all optimal algorithms for anticlustering 
#' that are available in anticlust.
#' 
#'@param x The data input. Can be one of two structures: (1) A feature
#'     matrix where rows correspond to elements and columns correspond
#'     to variables (a single numeric variable can be passed as a
#'     vector). (2) An N x N matrix dissimilarity matrix; can be an
#'     object of class \code{dist} (e.g., returned by
#'     \code{\link{dist}} or \code{\link{as.dist}}) or a \code{matrix}
#'     where the entries of the upper and lower triangular matrix
#'     represent pairwise dissimilarities.
#' @param K How many anticlusters should be created or alternatively:
#'     (a) A vector describing the size of each group (the latter
#'     currently only works for \code{objective = "dispersion")}.
#' @param objective The anticlustering objective, can be "diversity",
#'     "variance", "kplus" or "dispersion".
#' @param solver Optional. The solver used to obtain the optimal
#'     method.  Currently supports "glpk", "symphony", "lpSolve" and "gurobi". 
#'     See details.
#' @param time_limit Time limit in seconds, given to the solver.
#'    Default is there is no time limit.
#'
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and \code{K}) to each input element.
#' 
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @details 
#' 
#' This is a wrapper for all optimal methods supported in anticlust
#' (currently and in the future).  As compared to
#' \code{\link{anticlustering}}, it allows to specify the solver to
#' obtain an optimal solution and it can be used to obtain optimal
#' solutions for all supported anticlustering objectives (variance,
#' diversity, k-plus, dispersion). For the objectives "variance",
#' "diversity" and "kplus", the optimal ILP method in Papenberg and
#' Klau (2021) is used, which maximizes the sum of all pairwise
#' intra-cluster distances (given user specified number of clusters,
#' for equal-sized clusters).  To employ k-means anticlustering
#' (i.e. set \code{objective = "variance"}), the squared Euclidean
#' distance is used. For k-plus anticlustering, the squared Euclidean
#' distance based on the extended k-plus data matrix is used (see
#' \code{\link{kplus_moment_variables}}).  For the diversity (and the
#' dispersion), the Euclidean distance is used by default, but any
#' user-defined dissimilarity matrix is possible.
#' 
#' The dispersion is solved optimal using the approach described in
#' \code{\link{optimal_dispersion}}.
#' 
#' The optimal methods make use of "solvers" that actually implement
#' the algorithm for finding optimal solutions. The package anticlust
#' supports three solvers:
#' 
#' \itemize{
#'   \item{The default solver lpSolve (<https://sourceforge.net/projects/lpsolve/>).}
#'   \item{GNU linear programming kit (<http://www.gnu.org/software/glpk/>), 
#'   available from the package Rglpk and requested using \code{solver = "glpk"}.
#'   The R package Rglpk has to be installed manually if this solver should be used.}
#'   \item{The Symphony solver (<https://github.com/coin-or/SYMPHONY>),
#'   available from the package Rsymphony and requested using \code{solver = "symphony"}.
#'   (The package Rsymphony has to be installed manually if this solver should be used).}
#'   \item{The commercial gurobi solver, see https://www.gurobi.com/downloads/.}
#' }
#' 
#' For the maximum dispersion problem, it seems that the Symphony
#' solver is fastest, while the lpSolve solver seems to be good for
#' maximum diversity. However, note that in general the dispersion can
#' be solved optimally for much larger data sets than the diversity.
#' 
#' If a \code{time_limit} is set and the function cannot find in the optimal
#' objective in the given time, it will throw an error.
#' 
#' @export
#' 
#' @examples 
#' 
#' data <- matrix(rnorm(24), ncol = 2)
#' 
#' # These calls are equivalent for k-means anticlustering:
#' optimal_anticlustering(data, K = 2, objective = "variance")
#' optimal_anticlustering(dist(data)^2, K = 2, objective = "diversity")
#' 
#' # These calls are equivalent for k-plus anticlustering:
#' optimal_anticlustering(data, K = 2, objective = "kplus")
#' optimal_anticlustering(dist(kplus_moment_variables(data, 2))^2, K = 2, objective = "diversity")
#' 
optimal_anticlustering <- function(x, K, objective, solver = NULL, time_limit = NULL) {
  
  validate_input_optimal_anticlustering(x, K, objective, solver, time_limit)
  input_validation_anticlustering(
    x, K, objective = objective, method = "exchange", # 'method' is actually a lie to make it work
    preclustering = FALSE, categories = NULL, 
    repetitions = NULL, 
    standardize = FALSE
  )
  
  if (!argument_exists(solver)) {
    solver <- find_ilp_solver()
  }
  
  if (objective == "kplus") {
    x <- kplus_moment_variables(x, 2)
  }
  
  x <- convert_to_distances(x)
  if (objective %in% c("kplus", "variance")) {
    x <- x^2
    objective <- "diversity"
  }

  if (objective == "diversity") {
    ilp <- anticlustering_ilp(x, K)
    solution <- solve_ilp(ilp, solver = solver, time_limit = time_limit)
    if (solution$status != 0) {
      stop("Could not find the optimal objective in the given time limit.")
    }
    return(ilp_to_groups(solution, nrow(x)))
  } else {
    return(optimal_dispersion(x, K, solver, time_limit = time_limit)$groups)
  }
  
}
  
validate_input_optimal_anticlustering <- function(x, K, objective, solver, time_limit) {
  # x
  validate_data_matrix(x)
  # K
  validate_input(K, "K", objmode = "numeric", must_be_integer = TRUE, not_na = TRUE, not_function = TRUE)
  N <- nrow(as.matrix(x))
  if (!(length(K) == 1 || sum(K) == N)) {
    stop("Argument `K` is misspecified.")
  } 
  if (objective != "dispersion") {
    if (!all(K == K[1]) || length(K) == 1 && N%%K != 0) {
      stop("The selected method can only be used for equal-sized groups")
    }
  }
  
  # Objective
  validate_input(
    objective, "objective", objmode = "character", len = 1,
    input_set = c("variance", "diversity", "kplus", "dispersion"), 
    not_na = TRUE, 
    not_function = TRUE
  )

  # Solver
  if (argument_exists(solver)) {
    validate_input(solver, "solver", objmode = "character", len = 1,
                  input_set = c("glpk", "symphony", "lpSolve", "gurobi"), not_na = TRUE, not_function = TRUE)
    if (solver == "glpk") {
      if (!requireNamespace("Rglpk", quietly = TRUE)) {
        stop("The package Rglpk must be installed to use `solver = glpk`.\n", 
             "Type install.packages('Rglpk') in the R console to install it.")
      }
    } else if (solver == "symphony") {
      if (!requireNamespace("Rsymphony", quietly = TRUE)) {
        stop("The package Rsymphony must be installed to use `solver = symphony`.\n", 
             "Type install.packages('Rsymphony') in the R console to install it.")
      }
    } else if (solver == "gurobi") {
      if (!requireNamespace("gurobi", quietly = TRUE)) {
        stop("The package gurobi must be installed to use `solver = gurobi`.")
      }
    }
  }
  # time_limit 
  if (argument_exists(time_limit)) {
    validate_input(
      time_limit, 
      "time_limit",
      objmode = "numeric",
      len = 1, not_na = TRUE, 
      not_function = TRUE, 
      greater_than = 0,
      must_be_integer = TRUE
    )
  }
}
