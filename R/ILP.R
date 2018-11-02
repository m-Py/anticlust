
#' Construct the ILP represenation of an item assignment problem

#' @param distances An distance object or matrix representing the
#'   distances between items
#' @param p The number of groups to be created
#' @param solver A string identifing the solver to be used ("glpk",
#'   "gurobi", or "cplex")
#'
#' @return A list representing the ILP formulation of the instance
#'
#' @details To use this function, a linear programming solver must
#'  be installed and usable from R. The open source GNU linear
#'  programming kit (called from the package `Rglpk`) or one of the
#'  commercial solvers gurobi (called from the package `gurobi`) or
#'  IBM CPLEX (called from the package `Rcplex`) can be used. A license
#'  is needed for the commercial solvers. One of the interface packages
#'  must be installed.
#'
#' @importFrom Matrix Matrix
#'
#' @export
#'
#' @references
#'
#' T. Bulhões, G. F. de Sousa Filho, A. Subramanian, and F. C. Lucídio
#' dos Anjos, “Branch-and-cut approaches for p-cluster editing,”
#' Discrete Applied Mathematics, vol. 219, pp. 51–64, 2017.
#'
#' T. Bulhões, A. Subramanian, G. F. Sousa Filho, and F. C. Lucídio dos
#' Anjos, “Branch-and-price for p-cluster editing,” Computational
#' Optimization and Applications, vol. 67, no. 2, pp. 293–316, 2017.
#'
#' M. Grötschel and Y. Wakabayashi, “A cutting plane algorithm for a
#' clustering problem,” Mathematical Programming, vol. 45, nos. 1-3, pp.
#' 59–96, 1989.

item_assign_ilp <- function(distances, p, solver = "glpk") {

  ## identify solver because they use different identifiers for
  ## equality:
  if (solver == "glpk") {
    equal_sign <- "=="
    lower_sign <- "<="
    greater_sign <- ">="
  } else if (solver == "gurobi") {
    equal_sign <- "="
    lower_sign <- "<="
    greater_sign <- ">="
  } else if (solver == "cplex") {
    equal_sign <- "E"
    lower_sign <- "L"
    greater_sign <- "G"
  } else {
    stop("solver must be 'cplex', 'glpk', or 'gurobi'")
  }

  ## Problem: I have matrix of costs but need vector for ILP
  ## formulation.
  costs_m <- as.matrix(distances)
  n_items <-  nrow(costs_m)
  group_size = n_items / p
  ## Make vector of costs in data.frame (makes each cost identifiable)
  costs <- expand.grid(1:ncol(costs_m), 1:nrow(costs_m))
  colnames(costs) <- c("i", "j")
  costs$costs <- c(costs_m)
  ## TODO: add "_" to the names, i.e. xi_j, to avoid ambiguity when we
  ## have more than 100 items (which cannot be solved exactly
  ## probably, but we should have the option)
  costs$pair <- paste0("x", paste0(costs$i, costs$j))
  ## remove all cases where j >= i, i.e. remove redundant or self distances
  costs <- subset(costs, i < j)
  rownames(costs) <- NULL
  ## `costs` now contains the "vectorized" distances.

  ## How many decision variables do I have? Edges + coding of cluster
  ## leadership for each item.
  n_vars <- nrow(costs) + n_items

  ## Compute the number of constraints:
  ## The number of triangular constraints:
  n_tris <- choose(n_items, 3) * 3
  ## The number of constraints for constraint (5): The constraint
  ## numbers are taken from Bulhoes et al. 2017 ("Branch-and-price")
  n_c5 <- 1
  ## The number of constraints for constraint (6):
  n_c6 <- nrow(costs) ## for each decision var x_ij one constraint
  ## The number of constraints for constraint (7):
  n_c7 <- n_items - 1 ## for each but the first item one constraint
  ## The number of constraints for constraint (8):
  n_c8 <- 1
  ## Number of constraints enforcing the size of the groups:
  n_c9 <- n_items ## c9 not part of Bulhoes et al.'s formulation
  ## Total number of constraints:
  n_constraints <- n_tris + n_c5 + n_c6 + n_c7 + n_c8 + n_c9

  ## Start constructing the matrix representing the left-hand side of
  ## the constraints. Each column is a decision variable, each row is
  ## a constraint.

  constraints <- matrix(0, ncol = n_vars, nrow = n_constraints)
  colnames(constraints) <- c(costs$pair, paste0("y", 1:n_items))
  ## Use row and column names to identify the decision variables and
  ## the constraints. row names: "tc" are triangular constraints.
    rownames(constraints) <- c(paste0("tc", 1:n_tris), "c5",
                               paste0("c6_", 1:n_c6),
                               paste0("c7_", 1:n_c7),
                               "c8", paste0("c9_", 1:n_c9))

  ## (1) Triangular constraints

  counter <- 1
  for (i in 1:n_items) {
    for (j in 2:n_items) {
      for (k in 3:n_items) {
        ## ensure that only legal constraints are inserted:
        if (!(i < j) | !(j < k)) next
        ## offset for addressing the data.frame:
        offset <- (counter - 1) * 3
        ## triangular constraint 1
        constraints[offset + 1, paste0("x", i, j)] <- -1
        constraints[offset + 1, paste0("x", i, k)] <- 1
        constraints[offset + 1, paste0("x", j, k)] <- 1
        constraints[offset + 1, paste0("y", k)] <- 1
        ## triangular constraint 2
        constraints[offset + 2, paste0("x", i, j)] <- 1
        constraints[offset + 2, paste0("x", i, k)] <- -1
        constraints[offset + 2, paste0("x", j, k)] <- 1
        ## triangular constraint 3
        constraints[offset + 3, paste0("x", i, j)] <- 1
        constraints[offset + 3, paste0("x", i, k)] <- 1
        constraints[offset + 3, paste0("x", j, k)] <- -1
        ## increase counter
        counter <- counter + 1
      }
    }
  }
  constraints <- insert_group_contraints(constraints, n_items, i, j, group_size)


  ## Make the to-be-returned constraint matrix take less storage
  ## as a sparse matrix: (TODO: make it sparse from the beginning)
  constraints <- Matrix::Matrix(constraints, sparse = TRUE)

  ## (7) Insert the direction of the constraints:
  equalities = c(rep(lower_sign, n_tris),
                 rep(equal_sign, n_c5),
                 rep(lower_sign, n_c6),
                 rep(greater_sign, n_c7),
                 rep(equal_sign, n_c8),
                 rep(equal_sign, n_c9))

    ## (8) right hand side of the ILP <- many ones
    rhs <- c(rep(1, nrow(constraints) - 1 - n_c9),
             p, rep(group_size - 1, n_c9)) #  p = number of clusters


  ## (9) Construct objective function. Add values for the cluster leader
  ## decision variables to the objective functions. These must be 0
  ## because they are not part of the objective function, but they are
  ## needed for the standard form of the ILP
  obj_function <- c(costs$costs, rep(0, n_items))

  ## give names to all objects for inspection of the matrix
  names(rhs) <- rownames(constraints)
  names(equalities) <- rownames(constraints)
  names(obj_function) <- colnames(constraints)

  ## return instance
  instance              <- list()
  instance$n_groups     <- p
  instance$group_size   <- group_size
  instance$distances    <- distances
  instance$costs        <- costs
  instance$constraints  <- constraints
  instance$equalities   <- equalities
  instance$rhs          <- rhs
  instance$obj_function <- obj_function

  return(instance)
}

## This function inserts constraints concerned with cluster number and cluster
## size
insert_group_contraints <- function(constraints, n_items, i, j, group_size) {
  ## (2) Constraint (5): make first element leader of cluster
  constraints["c5", "y1"] <- 1

  ## (3) Constraint (6): restrictions on cluster leadership
  ## y_j <= 1 - x_ij <=> y_j + x_ij <= 1

  counter <- 1
  for (i in 1:n_items) {
    for (j in 2:n_items) {
      if (i >= j) next
      constraints[paste0("c6_", counter), paste0("x", i, j)] <- 1
      constraints[paste0("c6_", counter), paste0("y", j)] <- 1
      counter <- counter + 1
    }
  }

  ## (4) Constraint (7): more restrictions on cluster leadership
  ## y_j >= 1 - sum(x_ij) <=> y_j + sum(x_ij) >= 1
  counter <- 1
  for (j in 2:n_items) {
    constraints[paste0("c7_", counter), paste0("y", j)] <- 1
    for (i in 1:n_items) {
      if (i >= j) next
      constraints[paste0("c7_", counter), paste0("x", i, j)] <- 1
    }
    counter <- counter + 1
  }

  ## (5) Constraint (8): Specify number of clusters
  constraints["c8", paste0("y", 1:n_items)] <- 1

  ## (6) Add constraint forcing the clusters to be of equal size (I call
  ## it constraint (9) but it is not taken from Bulhoes 2017)

  for (i in 1:n_items) { ## i: current item that may be a cluster leader
    for (j in 1:n_items) {
      ## mark all edges that i is part of (i may be lower or higher index!)
      if (i < j)
        constraints[paste0("c9_", i), paste0("x", i, j)] <- 1
      if (j < i)
        constraints[paste0("c9_", i), paste0("x", j, i)] <- 1
    }
  }
  return(constraints)
}
