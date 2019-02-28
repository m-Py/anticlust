
#' Solve anticlustering using the distance criterion
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to items and columns correspond to features.
#' @param n_anticlusters How many anticlusters should be created.
#' @param solver A string identifying the solver to be used ("glpk",
#'     "gurobi", or "cplex")
#' @param standardize Boolean - should the feature values be
#'     standardized before groups are created?
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created.
#'
#' @return A vector representing the anticluster affiliation of
#'     elements.
#'
#' @noRd

exact_anticlustering <- function(features, n_anticlusters, solver, preclustering) {

  n_items <- nrow(features)
  distances <- dist(features)

  if (preclustering == TRUE) {
    ilp <- anticlustering_ilp(distances, n_items / n_anticlusters,
                              solver = solver)
    solution <- solve_ilp(ilp, solver, "min")
    assignment <- ilp_to_groups(ilp, solution)
    ## Fix distances - ensures that the most similar items are assigned
    ## to different groups
    distances <- edit_distances(distances, assignment)
    ## Edit ILP - objective function and group sizes
    ilp$obj_function <- vectorize_weights(distances)$costs
    ilp$rhs <- c(rep(1, choose(n_items, 3) * 3),
                 rep((n_items / n_anticlusters) - 1, n_items))
    ## Solve edited ILP
    solution <- solve_ilp(ilp, solver)
    assignment <- ilp_to_groups(ilp, solution)
    return(assignment)
  }

  ## Here the ILP is created without adjusting distances; i.e., true
  ## exact anticlustering
  ilp <- anticlustering_ilp(distances, n_anticlusters, solver = solver)
  solution <- solve_ilp(ilp, solver)
  ilp_to_groups(ilp, solution)
}


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
equal_sized_cluster_editing <- function(features, n_clusters, solver,
                                        standardize = FALSE) {

  if (standardize) {
    features <- scale(features)
  }
  distances <- dist(features)
  ilp <- anticlustering_ilp(distances, n_clusters, solver = solver)
  solution <- solve_ilp(ilp, solver, "min")
  assignment <- ilp_to_groups(ilp, solution)
  return(assignment)
}


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
                        control = list(round = 1, trace = 0)) # important to obtain integers
    ret_list$x <- ilp_solution$xopt
    ret_list$obj <- ilp_solution$obj
  }
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  return(ret_list)
}
