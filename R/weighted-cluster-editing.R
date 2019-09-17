
#' Solve weighted cluster editing exactly
#'
#' @param weights A N x N weight matrix. Larger values indicate stronger
#'     agreement (this is unlike in the \code{\link{anticlustering}} function)
#'
#' @return The clusters
#'
#' @examples
#' subs <- sample(150, size = 60)
#' data <- iris[subs, -5]
#' distances <- dist(data)
#' # Define agreement as being close enough to each other
#' agreements <- ifelse(as.matrix(distances) < 2.5, 1, -1)
#' clusters <- wce(agreements)
#' plot_clusters(data[, 1:2], clusters)
#'
#' @export
#'

wce <- function(weights) {
  solver <- solver_available()
  weights <- as.matrix(weights)
  ilp <- anticlustering_ilp(weights, K = 0, solver, FALSE) # k is irrelevant
  solution <- solve_ilp(ilp, solver, "max")
  ilp_to_groups(ilp, solution)
}
