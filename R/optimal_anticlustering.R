#' Optimal ("exact") algorithms for anticlustering
#' 
#' Wrapper function that gives access to all optimal algorithms for anticlustering 
#' that are available in anticlust.
#' 
#'@param x The data input. Can be one of two structures: (1) A
#'     feature matrix where rows correspond to elements and columns
#'     correspond to variables (a single numeric variable can be
#'     passed as a vector). (2) An N x N matrix dissimilarity matrix;
#'     can be an object of class \code{dist} (e.g., returned by
#'     \code{\link{dist}} or \code{\link{as.dist}}) or a \code{matrix}
#'     where the entries of the upper and lower triangular matrix
#'     represent pairwise dissimilarities. 
#' @param K How many anticlusters should be created or alternatively:
#'     (a) A vector describing the size of each group (the latter currently 
#'     only works for \code{objective = "dispersion")}.
#' @param objective The anticlustering objective, can be "diversity", "variance", 
#'     "kplus" or "dispersion".
#' @param solver Optional. The solver used to obtain the optimal method. 
#'     Currently supports "glpk" and "symphony". See details.
#'     
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and \code{K}) to each input element.
#' 
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @details 
#' 
#' This is a wrapper for all optimal methods supported in anticlust (currently and in the future). 
#' As compared to \code{\link{anticlustering}}, it allows to specify the solver to obtain an optimal
#' solution and it can be used to obtain optimal solutions for all supported
#' anticlustering objectives (variance, diversity, k-plus, dispersion). For 
#' the objectives "variance", "diversity" and "kplus", the optimal ILP method 
#' in Papenberg and Klau (2021) is used, which maximizes the sum of all pairwise 
#' intra-cluster distances (given user specified number of clusters, for equal-sized clusters).
#' To employ k-means anticlustering (i.e. set \code{objective = "variance"}), the
#' squared Euclidean distance is used. For k-plus anticlustering, the squared Euclidean distance
#' based on the extended k-plus data matrix is used (see \code{\link{kplus_moment_variables}}).
#' For the diversity (and the dispersion), the Euclidean distance is used by default, 
#' but any user-defined dissimilarity matrix is possible.
#' 
#' The dispersion is solved optimal using the approach described in \code{\link{optimal_dispersion}}.
#' 
#' The optimal methods either require the R package \code{Rglpk} and the GNU linear programming kit
#' (<http://www.gnu.org/software/glpk/>), or the R package
#' \code{Rsymphony} and the COIN-OR SYMPHONY solver libraries
#' (<https://github.com/coin-or/SYMPHONY>). If the argument \code{solver} is not 
#' specified by the user, the function will try to find the GLPK or SYMPHONY 
#' solver and throw an error if none is available. It will select the 
#' GLPK solver if both are available because some rare instances have been observed where
#' the SYMPHONY solver crashes on Mac computers. I would still try out the 
#' SYMPHONY solver to see if the unlikely crash occurs. However, this has to be 
#' set by the user (at least if both solver packages Rsymphony and Rglpk are available on the system).
#' 
#' @export
#' 
#' @examples 
#' 
#' # data <- matrix(rnorm(24), ncol = 2)
#' 
#' # These calls are equivalent for k-means anticlustering:
#' # optimal_anticlustering(data, K = 2, objective = "variance")
#' # optimal_anticlustering(dist(data)^2, K = 2, objective = "diversity")
#' 
#' # These calls are equivalent for k-plus anticlustering:
#' # optimal_anticlustering(data, K = 2, objective = "kplus")
#' # optimal_anticlustering(dist(kplus_moment_variables(data, 2))^2, K = 2, objective = "diversity")
#' 
optimal_anticlustering <- function(x, K, objective, solver = NULL) {
  
  validate_input_optimal_anticlustering(x, K, objective, solver)
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
    solution <- solve_ilp_diversity(ilp, solver = solver)
    return(ilp_to_groups(solution, nrow(x)))
  } else {
    return(optimal_dispersion(x, K, solver)$groups)
  }
  
}
  
validate_input_optimal_anticlustering <- function(x, K, objective, solver) {
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
                   input_set = c("glpk", "symphony"), not_na = TRUE, not_function = TRUE)
  }
  
}
