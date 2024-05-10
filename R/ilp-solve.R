
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

solve_ilp_diversity <- function(ilp, objective = "max", solver = NULL) {
  if (objective == "max") {
    max <- TRUE
  } else {
    max <- FALSE
  }
  
  if (is.null(solver)) {
    solver <- find_ilp_solver()
  }
  
  solver_function <- ifelse(solver == "symphony", Rsymphony::Rsymphony_solve_LP, Rglpk::Rglpk_solve_LP)
  name_opt <- ifelse(solver == "symphony", "objval", "optimum")
  
  ilp_solution <- solver_function(
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
  ret_list$obj <- ilp_solution[[name_opt]]
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}

# Function to find a solver package
find_ilp_solver <- function() {
  if (requireNamespace("Rglpk", quietly = TRUE)) {
    return("glpk")
  }
  else if (requireNamespace("Rsymphony", quietly = TRUE)) {
    return("symphony")
  }
  check_if_solver_is_available() # throws an error here!
}
