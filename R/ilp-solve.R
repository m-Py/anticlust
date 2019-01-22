
#' Solve an instance of item assignment
#'
#' @param items A data.frame of item features. Rows must correspond to
#'     items and columns to features.
#' @param n_groups How many groups are to be created.
#' @param solver A string identifying the solver to be used ("glpk",
#'   "gurobi", or "cplex")
#' @param standardize Boolean - should the feature values be
#'     standardized before groups are created? Defaults to FALSE.
#' @param heuristic Set the level of "heuristicism" by setting a numeric
#'   value of 0 to 3. Set to 0 to obtain the exact solution for the item
#'   assignment instance. Levels 1 to 3 will in a first step identify
#'   very similar items and then ensure that these will be assigned to
#'   different groups. Level 1 does this preclustering using exact
#'   cluster editing; Level 2 does this preclustering using a heuristic
#'   based on k-means clustering. On Level 3, preclustering is also done
#'   using the heuristic based on k-means clustering, and the final
#'   assignment is no longer done using exact ILP item assignment, but
#'   instead using a repeated random assignment.
#'
#' @return A vector representing the group assignment.
#'
#' @details To use this function, a linear programming solver must
#'  be installed and usable from R. The open source GNU linear
#'  programming kit (called from the package `Rglpk`) or one of the
#'  commercial solvers gurobi (called from the package `gurobi`) or
#'  IBM CPLEX (called from the package `Rcplex`) can be used. A license
#'  is needed for the commercial solvers. One of the interface packages
#'  must be installed.
#'
#' @export
#'
item_assignment <- function(items, n_groups, solver, standardize = FALSE,
                            heuristic = 1) {

  ## some input handling
  if (!(heuristic %in% 0:3)) {
    stop("argument heuristic must have a value between 0 and 3")
  }

  if (standardize) {
    items <- scale(items)
  }

  n_items <- nrow(items)
  distances <- dist(items)

  if (heuristic == 1) {
    ## Preclustering
    ilp <- item_assign_ilp(distances, n_items / n_groups,
                           solver = solver)
    solution <- solve_ilp(ilp, solver, "min")
    assignment <- ilp_to_groups(ilp, solution)
    ## Fix distances - ensures that the most similar items are assigned
    ## to different groups
    distances <- edit_distances(distances, assignment)
    ## Edit ILP - objective function and group sizes
    ilp$obj_function <- vectorize_weights(distances)$costs
    ilp$rhs <- c(rep(1, choose(n_items, 3) * 3),
                 rep((n_items / n_groups) - 1, n_items))
    ## Solve edited ILP to solve heuristic item assignment
    solution <- solve_ilp(ilp, solver)
    assignment <- ilp_to_groups(ilp, solution)
    return(assignment)
  } else if (heuristic > 1) {
    assignment <- equal_sized_clustering(data.frame(items), n_items / n_groups)
    distances  <- edit_distances(distances, assignment)
  }

  if (heuristic == 3) {
    assignment <- heuristic_item_assignment(items, assignment)
    return(assignment)
  }

  ## The following creates the item assignment ILP formulation and
  ## solves it:
  ilp <- item_assign_ilp(distances, n_groups, solver = solver)
  solution <- solve_ilp(ilp, solver)
  assignment <- ilp_to_groups(ilp, solution)
  return(assignment)
}


#' Edit distances of neighbours
#'
#' Based on a preclustering, distances between items of the same cluster
#' are set to a large negative value. By considering a preclustering of
#' items, very similar items will not be assigned to the same group when
#' the fixed distance object is used to create the ILP formulation of
#' the item assignment instance.
#'
#' @param distances A distance object or matrix of
#'     between-item-distances.
#' @return A vector representing the group assignment; objects that are
#'     part of the same group as indicated by this vector are assigned a
#'     new distance.
#' @param value The value that is assigned to the fixed distances.
#'     Defaults to -1,000,000 currently.
#'
#' @return A distance object containing the fixed distances. For items
#'     that are part of the same cluster (as specified in the
#'     `data.frame` `assignment`), the between-item-distances are set to
#'     -1000000. This will have to be replaced by a theoretically-sound
#'     value; -1000000 is just a hack that will work in the present
#'     applications.
#'
#' @importFrom stats as.dist dist
#' @export

edit_distances <- function(distances, assignment, value = -1000000) {
  n_groups <- length(unique(assignment))
  distances <- as.matrix(distances)
  for (i in 1:n_groups) {
    items <- which(assignment == i) ## which items are in the i'th cluster
    ## two for-loops to tap all pairwise distances that may require
    ## fixing (a lot of unnecessary iterations, probably)
    for (j in 1:(length(items) - 1)) {
      for (t in 2:length(items)) {
        distances[items[j], items[t]] <- value
        distances[items[t], items[j]] <- value
      }
    }
  }
  return(as.dist(distances))
}


#' Solve exact equal-sized cluster editing
#'
#' @param items A data.frame of item features. Rows must correspond to
#'     items and columns to features.
#' @param n_groups How many clusters are to be created.
#' @param solver A string identifing the solver to be used ("Rglpk",
#'   "gurobi", or "cplex")
#' @param standardize Boolean - should the feature values be
#'     standardized before groups are created? Defaults to FALSE.
#'
#' @return A vector representing the clustering.
#'
#' @export
#'
equal_sized_cluster_editing <- function(items, n_groups, solver,
                                        standardize = FALSE) {
  if (standardize) {
    items <- scale(items)
  }
  distances <- dist(items)
  ilp <- item_assign_ilp(distances, n_groups, solver = solver)
  solution <- solve_ilp(ilp, solver, "min")
  assignment <- ilp_to_groups(ilp, solution)
  return(assignment)
}


#' Solve the ILP formulation of the item assignment problem
#'
#' Usually it will be advised to call the higher level function
#' `item_assignment` By setting the `objective` parameter to "min", this
#' function solves weighted cluster editing instead of item assignment.
#'
#' @param ilp An object representing the ILP formulation of the
#'     instance, returned by `item_assign_ilp`
#' @param solver A string identifing the solver to be used ("glpk",
#'   "gurobi", or "cplex")
#' @param objective A string identifying whether the objective function
#'     of the ILP should be maximized ("max") or minimized
#'     ("min"). Maximizing creates similar groups (i.e., solves item
#'     assignment), minimizing creates distinct clusters (i.e., solves
#'     cluster editing).
#'
#' @return A `list` having two entries: `x` is the vector of optimal
#'   coefficients for all decision variables. `obj` is the optimal
#'   objective value.
#'
#' @details To use this function, a linear programming solver must
#'  be installed and usable from R. The open source GNU linear
#'  programming kit (called from the package `Rglpk`) or one of the
#'  commercial solvers gurobi (called from the package `gurobi`) or
#'  IBM CPLEX (called from the package `Rcplex`) can be used. A license
#'  is needed for the commercial solvers. One of the interface packages
#'  must be installed.
#'
#' @export
#'

solve_ilp <- function(ilp, solver, objective = "max") {

  ret_list <- list() # return the optimal value and the variable assignment

  if (solver == "glpk") {
    max <- FALSE
    if (objective == "max")
      max <- TRUE
    start <- Sys.time()
    ilp_solution <- Rglpk::Rglpk_solve_LP(obj = ilp$obj_function,
                                   mat = ilp$constraints,
                                   dir = ilp$equalities,
                                   rhs = ilp$rhs,
                                   types = "B",
                                   max = max)
    ret_list$x <- ilp_solution$solution
    ret_list$obj <- ilp_solution$optimum
    end <- Sys.time()
    print(end - start)
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
  } else if (solver == "cplex") {
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
