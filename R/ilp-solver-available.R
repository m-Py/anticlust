
#' Check if a solver package can be used
#'
#' Has no parameters, checks in the installed packages from the user
#'
#' @return \code{FALSE} If no solver is available; A string identifying
#'   the solver if at least one is available ("Rcplex", "gurobi", "Rglpk")
#'
#' @importFrom utils installed.packages
#'
#' @noRd
solver_available <- function() {
  solvers <- c("Rcplex", "gurobi", "Rglpk")
  pcks <- rownames(installed.packages())
  solvers_available <- solvers %in% pcks
  if (sum(solvers_available) == 0) # no solver available
    return(FALSE)
  return(solvers[solvers_available][1]) # pick only one solver
}

