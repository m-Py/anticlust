
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

  ## Initialize some constant variables:
  equality_signs <- equality_identifiers(solver)
  n_items <- nrow(as.matrix(distances))
  group_size = n_items / p
  costs <- vectorize_weights(distances)
  n_vars <- nrow(costs) + n_items # number of decision variables

  ## Specify the number of constraints:
  n_tris <- choose(n_items, 3) * 3 # tri-angular constraints
  n_c5 <- 1 # (constraint 5 in Bulhoes et al. 2017 ("Branch-and-price"))
  n_c6 <- nrow(costs) # numbering continues: 6, 7, ...
  n_c7 <- n_items - 1
  n_c8 <- 1
  n_c9 <- n_items ## constraint 9 not part of Bulhoes et al.'s formulation
  n_constraints <- n_tris + n_c5 + n_c6 + n_c7 + n_c8 + n_c9

  ## Construct ILP constraint matrix
  constraints <- matrix(0, ncol = n_vars, nrow = n_constraints)
  colnames(constraints) <- c(costs$pair, paste0("y", 1:n_items))
  rownames(constraints) <- c(paste0("tc", 1:n_tris), "c5",
                               paste0("c6_", 1:n_c6),
                               paste0("c7_", 1:n_c7),
                               "c8", paste0("c9_", 1:n_c9))
  constraints <- triangular_constraints(constraints, n_items)
  constraints <- group_contraints(constraints, n_items, group_size)
  constraints <- Matrix::Matrix(constraints, sparse = TRUE) ## TODO: init as sparse matrix


  ## Directions of the constraints:
  equalities <- c(rep(equality_signs$l, n_tris),
                  rep(equality_signs$e, n_c5),
                  rep(equality_signs$l, n_c6),
                  rep(equality_signs$g, n_c7),
                  rep(equality_signs$e, n_c8),
                  rep(equality_signs$e, n_c9))

  # Right-hand-side of ILP
  rhs <- c(rep(1, nrow(constraints) - 1 - n_c9),
           p, rep(group_size - 1, n_c9)) #  p = number of clusters

  # Objective function of the ILP
  obj_function <- c(costs$costs, rep(0, n_items))

  ## Give names to all objects for inspection purposes
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

#' Based on the solver, return identifiers for equality relationships
#' @param solver A string identifing the solver to be used ("Rglpk",
#'   "gurobi", or "cplex")
#'
#' @return A list of three elements containing strings representing
#'   equality (e), lower (l), and greater (g) relationships
#'
equality_identifiers <- function(solver) {
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
  return(list(e = equal_sign, l = lower_sign, g = greater_sign))
}


#' Function that returns a data.frame that contains vectorized weights
#'
#' @param distances A distance matrix
#' @return A data.frame having the following columns:
#'    `costs` - the actual weights in vectorized form
#'    `i` the first index of the item pair that is connected
#'    `j` the second index of the item pair that is connected
#'    `pair` A string of form x_ij identifying the item pair
vectorize_weights <- function(distances) {
  ## Problem: I have matrix of costs but need vector for ILP
  ## formulation.
  costs_m <- as.matrix(distances)
  ## Make vector of costs in data.frame (makes each cost identifiable)
  costs <- expand.grid(1:ncol(costs_m), 1:nrow(costs_m))
  colnames(costs) <- c("i", "j")
  costs$costs <- c(costs_m)
  ## remove redundant or self distances:
  costs <- subset(costs, i < j)
  ## TODO: add "_" to the names, i.e. xi_j, to avoid ambiguity when we
  ## have more than 100 items (which cannot be solved exactly
  ## probably, but we should have the option)
  costs$pair <- paste0("x", paste0(costs$i, costs$j))
  rownames(costs) <- NULL
  ## `costs` now contains the vectorized distances
  return(costs)
}

# Insert the coefficients of the triangular constraints into the constraint matrix
vectorized_triangular <- function(n_items, pair_names) {
  row_indices <- vector(length = choose(n_items, 3) * 10)
  col_indices <- row_indices
  xes         <- row_indices
  ## (1) Triangular constraints
  counter <- 1
  for (i in 1:n_items) {
    for (j in 2:n_items) {
      for (k in 3:n_items) {
        ## ensure that only legal constraints are inserted:
        if (!(i < j) | !(j < k)) next
        ## Offset for addressing the constraints
        offset_row <- (counter - 1) * 3 # for saying which row the constraint is in
        offset <- (counter - 1) * 10 # for targetting the vector
        ## triangular constraint 1
        row_indices[offset + c(1, 5, 8)]  <- offset_row + 1:3
        row_indices[offset + c(2, 6, 9)]  <- offset_row + 1:3
        row_indices[offset + c(3, 7, 10)] <- offset_row + 1:3
        row_indices[offset + 4] <- offset_row + 1

        col_indices[offset + c(1, 5, 8)] <- which(paste0("x", i, j) == pair_names)
        col_indices[offset + c(2, 6, 9)] <- which(paste0("x", i, k) == pair_names)
        col_indices[offset + c(3, 7, 10)] <- which(paste0("x", j, k) == pair_names)
        print(paste0("y", k))
        col_indices[offset + 4] <- which(paste0("y", k) == pair_names)

        xes[offset + c(1, 5, 8)] <- c(-1, 1, 1)
        xes[offset + c(2, 6, 9)] <- c(1, -1, 1)
        xes[offset + c(3, 7, 10)] <- c(1, 1, -1)
        xes[offset + 4] <- 1
        counter <- counter + 1
      }
    }
  }
  return(list(i = row_indices, j = col_indices, x = xes))
}


# Insert the coefficients of the triangular constraints into the constraint matrix
triangular_constraints <- function(constraints, n_items) {
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
  return(constraints)
}

## TODO next: vectorize group constraints (maybe only use group size constraint for a first test)
## Inserts constraints specifying cluster number and cluster size
group_contraints <- function(constraints, n_items, group_size) {
  ## The numbers are taken from Bulhoes et al. 2017 ("Branch-and-price")
  ## Constraint (5): make first element leader of cluster
  constraints["c5", "y1"] <- 1
  ## Constraint (6): restrictions on cluster leadership
  counter <- 1
  for (i in 1:n_items) {
    for (j in 2:n_items) {
      if (i >= j) next
      constraints[paste0("c6_", counter), paste0("x", i, j)] <- 1
      constraints[paste0("c6_", counter), paste0("y", j)] <- 1
      counter <- counter + 1
    }
  }
  ## Constraint (7): more restrictions on cluster leadership
  counter <- 1
  for (j in 2:n_items) {
    constraints[paste0("c7_", counter), paste0("y", j)] <- 1
    for (i in 1:n_items) {
      if (i >= j) next
      constraints[paste0("c7_", counter), paste0("x", i, j)] <- 1
    }
    counter <- counter + 1
  }
  ## Constraint (8): Specify number of clusters
  constraints["c8", paste0("y", 1:n_items)] <- 1
  ## Constraint (9) enforces the clusters to be of equal size (I call
  ## it constraint (9) but it is not taken from Bulhoes 2017)
  for (i in 1:n_items) { ## i: current item that may be a cluster leader
    for (j in 1:n_items) {
      if (i < j)
        constraints[paste0("c9_", i), paste0("x", i, j)] <- 1
      if (j < i) # i may be lower or higher index!)
        constraints[paste0("c9_", i), paste0("x", j, i)] <- 1
    }
  }
  return(constraints)
}
