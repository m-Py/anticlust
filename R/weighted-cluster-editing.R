
#' Exact weighted cluster editing
#'
#' @param weights A N x N weight matrix. Larger values indicate stronger
#'     agreement (this is unlike in the \code{\link{anticlustering}} function)
#'
#' @return An integer vector representing the cluster affiliation of each data point
#' 
#'
#' @examples
#' \dontrun{
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
#' To obtain an optimal solution for weighted cluster editing, 
#' a linear programming solver must be installed
#' and usable from R. The \code{anticlust} package 
#' requires the open source GNU linear programming kit, which 
#' is called from the package \code{Rglpk}).
#'

wce <- function(weights) {
  validate_solver()
  validate_data_matrix(weights)
  if (!is_distance_matrix(weights)) {
    stop("The input via argument `weights` is not a weight matrix, ",
         "the upper and lower triangulars of your matrix differ.")
  }
  weights <- as.matrix(weights)
  ilp <- anticlustering_ilp(weights, K = 0, FALSE) # k is irrelevant
  solution <- solve_ilp(ilp, "max")
  ilp_to_groups(solution, nrow(weights))
}
