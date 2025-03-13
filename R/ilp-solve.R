
#' Solve the ILP formulation of anticluster editing
#'
#' @param ilp An object representing the ILP formulation of the
#'     instance, returned by \code{anticlustering_ilp}
#' @param objective A string identifying whether the objective function
#'     of the ILP should be maximized ("max") or minimized
#'     ("min"). Maximizing creates similar groups (i.e., solves
#'     anticlustering), minimizing creates distinct clusters (i.e.,
#'     solves weighted cluster editing).
#' @param time_limit time limit given to solver, in seconds
#'
#' @return A `list` with two entries: `x` is the vector of optimal
#'     coefficients for all decision variables. `obj` is the optimal
#'     objective value.
#'
#' @noRd

solve_ilp <- function(ilp, objective = "max", solver = NULL, time_limit = NULL) {

  if (is.null(solver)) {
    solver <- find_ilp_solver()
  }

  if (solver == "glpk") {
    return(solve_ilp_glpk(ilp, objective, time_limit))
  } else if (solver == "symphony") {
    return(solve_ilp_symphony(ilp, objective, time_limit))
  } else if (solver == "lpSolve") {
    return(solve_ilp_lpSolve(ilp, objective, time_limit))
  } else if (solver == "gurobi") {
    return(solve_ilp_gurobi(ilp, objective, time_limit))
  }
}

solve_ilp_glpk <- function(ilp, objective, time_limit) {
  ilp_solution <- Rglpk::Rglpk_solve_LP(
    obj = ilp$obj_function,
    mat = ilp$constraints,
    dir = ilp$equalities,
    rhs = ilp$rhs,
    types = "B",
    max = objective == "max",
    control = list(tm_limit = ifelse(is.null(time_limit), 0, time_limit * 1000))
  )
  
  # return the optimal value and the variable assignment
  ret_list <- list() 
  ret_list$x <- ilp_solution$solution
  ret_list$obj <- ilp_solution$optimum
  ret_list$status <- ilp_solution$status
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}

solve_ilp_symphony <- function(ilp, objective, time_limit) {
  
  ilp_solution <- Rsymphony::Rsymphony_solve_LP(
    obj = ilp$obj_function,
    mat = ilp$constraints,
    dir = ilp$equalities,
    rhs = ilp$rhs,
    types = "B",
    max = objective == "max",
    time_limit = ifelse(is.null(time_limit), -1, time_limit)
  )
  
  # return the optimal value and the variable assignment
  ret_list <- list() 
  ret_list$x <- ilp_solution$solution
  ret_list$obj <- ilp_solution$objval
  ret_list$status <- ilp_solution$status
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}

solve_ilp_lpSolve <- function(ilp, objective, time_limit) {
  ilp_solution <- lpSolve::lp(
    objective,
    ilp$obj_function,
    as.matrix(ilp$constraints),
    ilp$equalities,
    ilp$rhs,
    all.bin = TRUE,
    timeout = ifelse(is.null(time_limit), 0, time_limit)
  )
  # return the optimal value and the variable assignment
  ret_list <- list() 
  ret_list$x <- ilp_solution$solution
  ret_list$obj <- ilp_solution$objval
  ret_list$status <- ilp_solution$status
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}

solve_ilp_gurobi <- function(ilp, objective, time_limit) {
  
  ilp$equalities[ilp$equalities == "=="] <- "=" # gurobi uses 1 '='

  ## build model
  model <- list()
  model$A          <- ilp$constraints
  model$obj        <- ilp$obj_function
  model$modelsense <- objective
  model$rhs        <- ilp$rhs
  model$sense      <- ilp$equalities
  model$vtype     <- "B"
  
  ## solve
  if (argument_exists(time_limit)) {
    ilp_solution <- gurobi::gurobi(model, params = list(LogToConsole = 0, TimeLimit = time_limit, MIPGap = 0))
  } else {
    ilp_solution <- gurobi::gurobi(model, params = list(LogToConsole = 0, MIPGap = 0))
  }
  ret_list <- list() 
  ret_list$x <- ilp_solution$x
  ret_list$obj <- ilp_solution$objval
  ret_list$status <- ifelse(ilp_solution$status == "OPTIMAL", 0, 1)
  ## name the decision variables
  if (!is.null(ilp_solution$x)) { # gurobi does not return solution if not optimal
    names(ret_list$x) <- colnames(ilp$constraints)
  }
  ret_list
}

# Function to find a solver package
find_ilp_solver <- function() {
  if (requireNamespace("gurobi", quietly = TRUE)) {
    return("gurobi")
  } else if (requireNamespace("lpSolve", quietly = TRUE)) {
    return("lpSolve")
  } else if (requireNamespace("Rglpk", quietly = TRUE)) {
    return("glpk")
  } else if (requireNamespace("Rsymphony", quietly = TRUE)) {
    return("symphony")
  }
  check_if_solver_is_available() # throws an error here!
}
