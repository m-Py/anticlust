
#' Solve anticlustering using the distance criterion
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to items and columns correspond to features.
#' @param n_anticlusters How many anticlusters should be created.
#' @param solver A string identifying the solver to be used ("glpk",
#'     "gurobi", or "cplex")
#' @param standardize Boolean - should the feature values be
#'     standardized before groups are created? Defaults to TRUE
#' @param heuristic Set the level of "heuristicism" by setting a numeric
#'     value of 0 to 3. Set to 0 to obtain the exact solution for the
#'     item assignment instance. Levels 1 to 3 will in a first step
#'     identify very similar items and then ensure that these will be
#'     assigned to different groups. Level 1 does this preclustering
#'     using exact cluster editing; Level 2 does this preclustering
#'     using a heuristic based on k-means clustering. On Level 3,
#'     preclustering is also done using the heuristic based on k-means
#'     clustering, and the final assignment is no longer done using
#'     exact ILP item assignment, but instead using a repeated random
#'     assignment.
#'
#' @return A vector representing the anticluster affiliation of
#'     elements.
#'
#' @details This function includes an exact approach to anticlustering
#'     using the distance criterion. To use this functionality, a linear
#'     programming solver must be installed and usable from R. The open
#'     source GNU linear programming kit (called from the package
#'     `Rglpk`) or one of the commercial solvers gurobi (called from the
#'     package `gurobi`) or IBM CPLEX (called from the package `Rcplex`)
#'     can be used. A license is needed for the commercial solvers. One
#'     of the interface packages must be installed. To improve run time,
#'     some restrictions can be enforced (see argument heuristic); the
#'     solution is still found through an ILP but the solution need not
#'     be optimal, though it often still is.
#'
#'
#' @export
#'
distance_anticlustering <- function(features, n_anticlusters, solver,
                                    standardize = TRUE, heuristic = 1) {

  ## some input handling
  if (!(heuristic %in% 0:3)) {
    stop("argument heuristic must have a value between 0 and 3")
  }

  if (standardize) {
    features <- scale(features)
  }

  n_items <- nrow(features)
  distances <- dist(features)

  if (heuristic == 1) {
    ## Preclustering
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
    ## Solve edited ILP to solve heuristic item assignment
    solution <- solve_ilp(ilp, solver)
    assignment <- ilp_to_groups(ilp, solution)
    return(assignment)
  } else if (heuristic > 1) {
    assignment <- equal_sized_kmeans(data.frame(features), n_items / n_anticlusters)
    distances  <- edit_distances(distances, assignment)
  }

  if (heuristic == 3) {
    assignment <- heuristic_anticlustering(features, assignment)
    return(assignment)
  }

  ## The following creates the item assignment ILP formulation and
  ## solves it:
  ilp <- anticlustering_ilp(distances, n_anticlusters, solver = solver)
  solution <- solve_ilp(ilp, solver)
  assignment <- ilp_to_groups(ilp, solution)
  return(assignment)
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
#' @export
#'
equal_sized_cluster_editing <- function(features, n_clusters, solver,
                                        standardize = FALSE) {

  if (standardize) {
    features <- scale(features)
  }
  distances <- dist(features)
  ilp <- distance_anticlustering(distances, n_clusters, solver = solver)
  solution <- solve_ilp(ilp, solver, "min")
  assignment <- ilp_to_groups(ilp, solution)
  return(assignment)
}


#' Solve the ILP formulation of distance anticlustering
#'
#' Usually it will be advised to call the higher level function
#' `distance_anticlustering`. By setting the `objective` parameter to
#' "min", this function solves weighted cluster editing instead of
#' anticlustering.
#'
#' @param ilp An object representing the ILP formulation of the
#'     instance, returned by `item_assign_ilp`
#' @param solver A string identifing the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
#' @param objective A string identifying whether the objective function
#'     of the ILP should be maximized ("max") or minimized
#'     ("min"). Maximizing creates similar groups (i.e., solves item
#'     assignment), minimizing creates distinct clusters (i.e., solves
#'     cluster editing).
#'
#' @return A `list` having two entries: `x` is the vector of optimal
#'     coefficients for all decision variables. `obj` is the optimal
#'     objective value.
#'
#' @details To use this function, a linear programming solver must be
#'     installed and usable from R. The open source GNU linear
#'     programming kit (called from the package `Rglpk`) or one of the
#'     commercial solvers gurobi (called from the package `gurobi`) or
#'     IBM CPLEX (called from the package `Rcplex`) can be used. A
#'     license is needed for the commercial solvers. One of the
#'     interface packages must be installed.
#'
#' @export
#'

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
                        control = list(round = 1)) # important to obtain integers
    ret_list$x <- ilp_solution$xopt
    ret_list$obj <- ilp_solution$obj
  }
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  return(ret_list)
}