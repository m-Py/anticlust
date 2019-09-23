
#' Solve exact equal-sized cluster editing
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to items and columns correspond to features.
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
