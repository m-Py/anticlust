
#' Exact weighted cluster editing
#'
#' Optimally solves weighted cluster editing (also known as »correlation clustering« or
#' »clique partitioning problem«). 
#'
#' @param x A N x N similarity matrix. Larger values indicate stronger
#'     agreement / similarity between a pair of data points
#' @param solver Optional argument; if passed, has to be either "glpk" or
#'   "symphony". See details.
#' @param solve_sequentially Logical (default FALSE). Should the ILP be solved 
#' sequentially, meaning that triangular constraints are only added to the model
#' if they are violated, and then the ILP is solved again.
#'   
#' @return An integer vector representing the cluster affiliation of each data point
#' 
#'
#' @examples
#' \donttest{
#' features <- swiss
#' distances <- dist(scale(swiss))
#' hist(distances)
#' # Define agreement as being close enough to each other.
#' # By defining low agreement as -1 and high agreement as +1, we
#' # solve *unweighted* cluster editing
#' agreements <- ifelse(as.matrix(distances) < 3, 1, -1)
#' clusters <- wce(agreements)
#' plot(swiss, col = clusters, pch = 19)
#' clusters2 <- wce(agreements, solve_sequentially = TRUE)
#' clusters1 == clusters2
#' }
#' 
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @details 
#' 
#' Finds the clustering that maximizes the sum of pairwise similarities within clusters
#' according to the model by Grötschel and Wakabayashi (1989), using integer linear
#' programming (ILP). In the input some similarities should be negative (indicating dissimilarity) because 
#' otherwise the maximum sum of similarities is obtained by simply joining all elements 
#' within a single big cluster. If the argument \code{solver} is not specified, the function
#' will try to find the GLPK or SYMPHONY solver by itself (it prioritizes using SYMPHONY if both are
#' available).
#' 
#' @note
#' 
#' This function either requires the R package \code{Rglpk} and the GNU linear 
#' programming kit (<http://www.gnu.org/software/glpk/>) or the R package 
#' \code{Rsymphony} and the COIN-OR SYMPHONY solver libraries 
#' (<https://github.com/coin-or/SYMPHONY>).
#' 
#' @references
#'
#' 
#' Bansal, N., Blum, A., & Chawla, S. (2004). Correlation clustering. 
#' Machine Learning, 56, 89–113. 
#' 
#' Böcker, S., & Baumbach, J. (2013). Cluster editing. In Conference on 
#' Computability in Europe (pp. 33–44).
#' 
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59-96.
#' 
#' Wittkop, T., Emig, D., Lange, S., Rahmann, S., Albrecht, M., Morris, J. H., ..., Baumbach,
#' J. (2010). Partitioning biological data with transitivity clustering. Nature Methods, 7,
#' 419–420. 
#'

wce <- function(x, solver = NULL, solve_sequentially = FALSE) {

  if (argument_exists(solver)) {
    validate_input(solver, "solver", objmode = "character", len = 1,
                   input_set = c("glpk", "symphony"), not_na = TRUE, not_function = TRUE)
  } else {
    solver <- find_ilp_solver()
  }

  validate_data_matrix(x)
  if (!is_distance_matrix(x)) {
    stop("The input via argument `weights` is not a similarity matrix, ",
         "the upper and lower triangulars of your matrix differ.")
  }
  
  x <- as.matrix(x)
  
  if (solve_sequentially) {
    solution <- wce_solve_sequentially(x, solver)
  } else {
    ilp <- anticlustering_ilp(x, K = 0, FALSE) # k is irrelevant
    solution <- solve_ilp_diversity(ilp, "max", solver)
  }
  ilp_to_groups(solution, nrow(x))
}

# Solve clique partitioning problem by subsequently solving instances and
# adding only violated triangular constraints
wce_solve_sequentially <- function(x, solver = NULL) {
  if (is.null(solver)) {
    solver <- find_ilp_solver()
  }
  solver_function <- ifelse(
    solver == "symphony", 
    Rsymphony::Rsymphony_solve_LP, 
    Rglpk::Rglpk_solve_LP
  )
  name_opt <- ifelse(
    solver == "symphony", 
    "objval", 
    "optimum"
  )
  # Initialize some variables:
  N <- nrow(x)
  costs <- vectorize_weights(x)
  N_VARS <- length(costs$i)
  # initialize constraint matrix with 1 constraint (between elements 1, 2, and 3)
  constraints_init <- Matrix::sparseMatrix(
    i = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
    j = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
    x = c(-1, 1, 1, 1, -1, 1, 1, 1, -1),
    dims = c(3, N_VARS)
  )
  OBJECTIVE <- costs$costs
  names(OBJECTIVE) <- costs$pair
  constraintID <- 1:length(costs$pair)
  names(constraintID) <- costs$pair
  constraints <- constraints_init
  
  violated <- TRUE
  
  while (violated) {
    N_CONSTRAINTS <- nrow(constraints)
    ilp_solution <- solver_function(
      obj = OBJECTIVE,
      mat = constraints,
      dir = rep("<=", N_CONSTRAINTS),
      rhs = rep(1, N_CONSTRAINTS),
      types = "B",
      max = TRUE
    )
    
    SOLUTION <- ilp_solution$solution
    names(SOLUTION) <- costs$pair
    
    violations <- c()
    
    # check all transitivity constraints
    for (i in 1:(N-2)) {
      for (j in (i+1):(N-1)) {
        for (k in (j+1):N) {
          # investigate the connections between the triplet
          c1 <- SOLUTION[paste0("x", i, "_", j)] 
          c2 <- SOLUTION[paste0("x", i, "_", k)] 
          c3 <- SOLUTION[paste0("x", j, "_", k)]
          sum1 <- c1 + c2 - c3
          sum2 <- c1 - c2 + c3
          sum3 <- -c1 + c2 + c3
          # add indices of violated constraints
          if (any(c(sum1, sum2, sum3) > 1)) {
            ID1 <- constraintID[paste0("x", i, "_", j)] 
            ID2 <- constraintID[paste0("x", i, "_", k)] 
            ID3 <- constraintID[paste0("x", j, "_", k)]
            violations <- c(violations, ID1, ID2, ID3, ID1, ID2, ID3, ID1, ID2, ID3)
          }
        }
      }
    }
    violated <- length(violations) > 0
    if (violated == TRUE) {
      row_ids <- rep(1:(length(violations)/3), each = 3)
      xes <- rep_len(c(-1, 1, 1, 1, -1, 1, 1, 1, -1), length(violations))
      constraints <- rbind(
        constraints, 
        sparseMatrix(
          i = row_ids,
          j = violations,
          x = xes,
          dims = c(length(violations)/3, N_VARS)
        )
      )
    } else {
      # return the optimal value and the variable assignment
      ret_list <- list() 
      ret_list$x <- ilp_solution$solution
      ret_list$obj <- ilp_solution[[name_opt]]
      ## name the decision variables
      names(ret_list$x) <- costs$pair
      return(ret_list)
    }
  }
  stop("THIS SHOULD NOT HAVE HAPPENED, PLEASE CONTACT THE ANTICLUST MAINTAINER.")
}


