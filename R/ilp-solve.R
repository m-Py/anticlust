
#' Solve the ILP formulation of anticluster editing
#'
#' @param ilp An object representing the ILP formulation of the
#'     instance, returned by \code{anticlustering_ilp}
#' @param objective A string identifying whether the objective function
#'     of the ILP should be maximized ("max") or minimized
#'     ("min"). Maximizing creates similar groups (i.e., solves
#'     anticlustering), minimizing creates distinct clusters (i.e.,
#'     solves weighted cluster editing).
#'
#' @return A `list` with two entries: `x` is the vector of optimal
#'     coefficients for all decision variables. `obj` is the optimal
#'     objective value.
#'
#' @noRd

solve_ilp <- function(ilp, objective = "max") {
  if (objective == "max") {
    max <- TRUE
  } else {
    max <- FALSE
  }
  ilp_solution <- Rglpk::Rglpk_solve_LP(
    obj = ilp$obj_function,
    mat = ilp$constraints,
    dir = ilp$equalities,
    rhs = ilp$rhs,
    types = "B",
    max = max
  )
  # return the optimal value and the variable assignment
  ret_list <- list() 
  ret_list$x <- ilp_solution$solution
  ret_list$obj <- ilp_solution$optimum
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}
