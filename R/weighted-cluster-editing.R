
#' Exact weighted cluster editing
#'
#' Optimally solves weighted cluster editing (also known as »correlation clustering« or
#' »clique partitioning problem«). 
#'
#' @param x A N x N similarity matrix. Larger values indicate stronger
#'     agreement / similarity between a pair of data points
#' @param K Optional argument; if specified, denotes the number of (equal-sized) clusters.
#' @param solver Either "glpk" (default) or "symphony"
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
#' }
#' 
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @details 
#' 
#' Finds the clustering that maximizes the sum of pairwise similarities within clusters. 
#' In the input some similarities should be negative (indicating dissimilarity) because 
#' otherwise the maximum sum of similarities is obtained by simply joining all elements 
#' within a single big cluster.
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
#' Wittkop, T., Emig, D., Lange, S., Rahmann, S., Albrecht, M., Morris, J. H., . . . Baumbach,
#' J. (2010). Partitioning biological data with transitivity clustering. Nature Methods, 7,
#' 419–420. 
#'

wce <- function(x, K = NULL, solver = NULL) {

  if (argument_exists(K)) {
    validate_input(K, "K", objmode = "numeric", must_be_integer = TRUE, not_na = TRUE, not_function = TRUE)
    equal_sized_groups <- TRUE
  } else {
    equal_sized_groups <- FALSE
    K <- 0 # k is irrelevant for standard WCE
  }
  
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
  ilp <- anticlustering_ilp(x, K = K, equal_sized_groups)
  solution <- solve_ilp_diversity(ilp, "max", solver)
  ilp_to_groups(solution, nrow(x))
}
