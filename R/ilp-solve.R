
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
  ## build model
  model <- list()
  model$A          <- ilp$constraints
  model$obj        <- ilp$obj_function
  model$modelsense <- objective
  model$rhs        <- ilp$rhs
  model$sense      <- ilp$equalities
  model$vtypes     <- "B"
  ## solve
  ilp_solution <- gurobi::gurobi(model, params = list(LogToConsole = 0))
  ret_list <- list()
  ret_list$x <- ilp_solution$x
  ret_list$obj <- ilp_solution$objval
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}
