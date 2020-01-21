
#' Exact weighted cluster editing
#'
#' @param weights A N x N weight matrix. Larger values indicate stronger
#'     agreement (this is unlike in the \code{\link{anticlustering}} function)
#'
#' @return The clusters
#'
#' @examples
#' features <- swiss
#' distances <- dist(scale(swiss))
#' hist(distances)
#' # Define agreement as being close enough to each other.
#' # By defining low agreement as -1 and high agreement as +1, we
#' # solve *unweighted* cluster editing
#' agreements <- ifelse(as.matrix(distances) < 3, 1, -1)
#' clusters <- wce(agreements)
#' plot(swiss, col = clusters, pch = 19)
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'

wce <- function(weights) {
  if (!is_distance_matrix(weights)) {
    stop("The input via argument `weights` is not a distance matrix. ",
         "Maybe the upper and lower triangulars of your matrix differ.")
  }
  solver <- solver_available()
  weights <- as.matrix(weights)
  ilp <- anticlustering_ilp(weights, K = 0, solver, FALSE) # k is irrelevant
  solution <- solve_ilp(ilp, solver, "max")
  ilp_to_groups(solution, nrow(weights))
}
