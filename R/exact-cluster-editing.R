
#' Solve exact equal-sized cluster editing
#'
#' @param data A N x N distance matrix or N x M features table.
#' @param K How many clusters are to be created.
#' @param solver A string identifing the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
#'
#' @return A vector representing the clustering.
#'
#' @noRd
#'
balanced_cluster_editing <- function(data, K, solver) {
  if (!is_distance_matrix(data)) {
    data <- as.matrix(dist(data))
  }
  ilp <- anticlustering_ilp(data, K, solver = solver)
  solution <- solve_ilp(ilp, solver, "min")
  ilp_to_groups(solution, nrow(data))
}
