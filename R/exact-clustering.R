
#' Solve exact equal-sized cluster editing
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to items and columns correspond to features.
#' @param n_clusters How many clusters are to be created.
#' @param solver A string identifing the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
#' @param standardize Boolean, should the feature values be standardized
#'     before groups are created? Defaults to FALSE.
#'
#' @return A vector representing the clustering.
#'
#' @noRd
#'
balanced_cluster_editing <- function(distances, n_clusters, solver,
                                     standardize = FALSE) {
  ilp <- anticlustering_ilp(distances, n_clusters, solver = solver)
  solution <- solve_ilp(ilp, solver, "min")
  ilp_to_groups(ilp, solution)
}
