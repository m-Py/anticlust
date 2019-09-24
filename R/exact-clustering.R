
#' Solve exact equal-sized cluster editing
#'
#' @param distances A N x N distance matrix (a «matrix», not of type dist)
#' @param K How many clusters are to be created.
#' @param solver A string identifing the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
#'
#' @return A vector representing the clustering.
#'
#' @noRd
#'
balanced_cluster_editing <- function(distances, K, solver) {
  ilp <- anticlustering_ilp(distances, K, solver = solver)
  solution <- solve_ilp(ilp, solver, "min")
  ilp_to_groups(ilp, solution)
}
