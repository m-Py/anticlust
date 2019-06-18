
#' Solve the ILP formulation of distance anticlustering
#'
#'
#' @param ilp An object representing the ILP formulation of the
#'     instance, returned by \code{anticlustering_ilp}
#' @param solver A string identifing the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
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

solve_ilp <- function(ilp, solver, objective = "max") {

  ret_list <- list() # return the optimal value and the variable assignment

  if (solver == "Rglpk") {
    max <- FALSE
    if (objective == "max")
      max <- TRUE
    ilp_solution <- Rglpk::Rglpk_solve_LP(obj = ilp$obj_function,
                                          mat = ilp$constraints,
                                          dir = ilp$equalities,
                                          rhs = ilp$rhs,
                                          types = "B",
                                          max = max)
    ret_list$x <- ilp_solution$solution
    ret_list$obj <- ilp_solution$optimum
  } else if (solver == "gurobi") {
    ## build model
    model <- list()
    model$A          <- ilp$constraints
    model$obj        <- ilp$obj_function
    model$modelsense <- objective
    model$rhs        <- ilp$rhs
    model$sense      <- ilp$equalities
    model$vtypes     <- "B"
    ## solve
    ilp_solution <- gurobi::gurobi(model)
    ret_list$x <- ilp_solution$x
    ret_list$obj <- ilp_solution$objval
  } else if (solver == "Rcplex") {
    ilp_solution <- Rcplex::Rcplex(cvec = ilp$obj_function,
                                   Amat = ilp$constraints,
                                   bvec = ilp$rhs,
                                   objsense = objective,
                                   sense = ilp$equalities,
                                   vtype = "I",
                                   control = list(round = 1, trace = 0, epgap = 1e-11)) # important to obtain integers
    ret_list$x <- ilp_solution$xopt
    ret_list$obj <- ilp_solution$obj
  }
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  return(ret_list)
}
