
#' Solve exact equal-sized cluster editing
#'
#' @param data A N x N distance matrix or N x M features table.
#' @param K How many clusters are to be created.
#' @param solver "glpk" or "symphony"
#'
#' @return A vector representing the clustering.
#'
#' @noRd
#' 
#'
balanced_cluster_editing <- function(data, K, solver) {
  if (!is_distance_matrix(data)) {
    data <- as.matrix(dist(data))
  }
  ilp <- anticlustering_ilp(data, K)
  solution <- solve_ilp_diversity(ilp, "min", solver)
  ilp_to_groups(solution, nrow(data))
}
