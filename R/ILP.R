
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
#' @export
#'
#' @references
#'
#' M. Grötschel and Y. Wakabayashi, “A cutting plane algorithm for a
#' clustering problem,” Mathematical Programming, vol. 45, nos. 1-3, pp.
#' 59–96, 1989.

item_assign_ilp <- function(distances, p, solver = "glpk") {

  ## Initialize some constant variables:
  equality_signs <- equality_identifiers(solver)
  n_items        <- nrow(as.matrix(distances))
  group_size     <- n_items / p
  costs          <- vectorize_weights(distances)

  ## Specify the number of constraints:
  n_tris <- choose(n_items, 3) * 3 # triangular constraints
  n_c9 <- n_items # number of group constraints

  ## Construct ILP constraint matrix
  constraints <- sparse_constraints(n_items, costs$pair)
  colnames(constraints) <- costs$pair

  ## Directions of the constraints:
  equalities <- c(rep(equality_signs$l, n_tris),
                  rep(equality_signs$e, n_c9))

  # Right-hand-side of ILP
  rhs <- c(rep(1, n_tris), rep(group_size - 1, n_c9)) #  p = number of clusters

  # Objective function of the ILP
  obj_function <- costs$costs

  ## Give names to all objects for inspection purposes
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


#' Convert matrix of distances into vector of distances
#'
#' @param distances A distance matrix
#' @return A data.frame having the following columns:
#'    `costs` - the actual weights in vectorized form
#'    `i` the first index of the item pair that is connected
#'    `j` the second index of the item pair that is connected
#'    `pair` A string of form x_ij identifying the item pair
vectorize_weights <- function(distances) {
  ## Problem: I have matrix of costs but need vector for ILP
  costs_m <- as.matrix(distances)
  ## Make vector of costs in data.frame (makes each cost identifiable)
  costs <- expand.grid(1:ncol(costs_m), 1:nrow(costs_m))
  colnames(costs) <- c("i", "j")
  costs$costs <- c(costs_m)
  ## remove redundant or self distances:
  costs <- subset(costs, i < j)
  costs$pair <- paste0("x", paste0(costs$i, "_", costs$j))
  rownames(costs) <- NULL
  return(costs)
}

#' Construct a sparse matrix representing the ILP constraints
#' @param n_items How many items are there
#' @param pair_names A character vector of names representing the item pairs.
#'   Must have the form that is contained in the ILP data.frame costs$pair.
#' @return A sparse matrix representing the left-hand side of the ILP (A in Ax ~ b)
#'
#' @importFrom Matrix sparseMatrix
sparse_constraints <- function(n_items, pair_names) {
  tri <- vectorized_triangular(n_items, pair_names)
  gr  <- vectorized_group(n_items, pair_names)
  return(Matrix::sparseMatrix(c(tri$i, gr$i), c(tri$j, gr$j), x = c(tri$x, gr$x)))
}

#' Construct indices for a sparse matrix representation of triangular constraints
#' @param n_items How many items are there
#' @param pair_names A character vector of names representing the item pairs.
#'   Must have the form that is contained in the ILP data.frame costs$pair.
#' @return A list of indices to be used as input parameters of Matrix::sparseMatrix
#'
vectorized_triangular <- function(n_items, pair_names) {
  triangular_constraints <- choose(n_items, 3)
  coef_per_constraint <- 3
  # number of coefficients in constraint matrix:
  col_indices <- matrix(ncol = triangular_constraints, nrow = coef_per_constraint * 3)
  row_indices <- rep(1:(triangular_constraints*3), each = 3)
  xes <- rep(c(-1, 1, 1, 1, -1, 1, 1, 1, -1), triangular_constraints)
  ## Fill columns for constraints
  counter <- 1
  for (i in 1:n_items) {
    for (j in 2:n_items) {
      for (k in 3:n_items) {
        ## ensure that only legal constraints are inserted:
        if (!(i < j) | !(j < k)) next
        ## Construct indices
        pairs <- c(paste0("x", i, "_", j), paste0("x", i, "_", k), paste0("x", j, "_", k))
        indices <- match(pairs, pair_names)
        col_indices[, counter] <- rep(indices, 3)
        counter <- counter + 1
      }
    }
  }
  return(list(i = row_indices, j = c(col_indices), x = xes))
}

#' Construct indices for a sparse matrix representation of group constraints
#' @param n_items How many items does the instance have
#' @param pair_names A character vector of names representing the item pairs.
#'   Must have the form that is contained in the ILP data.frame costs$pair.
#' @return A list of indices to be used as input parameters of Matrix::sparseMatrix
vectorized_group <- function(n_items, pair_names) {
  coef_per_constraint <- (n_items - 1)
  group_coefficients <- coef_per_constraint * n_items
  row_indices <- rep((1:n_items) + (3 * choose(n_items, 3)), each = coef_per_constraint)
  col_indices <- matrix(ncol = n_items, nrow = coef_per_constraint)
  xes         <- rep(1, length = group_coefficients)

  for (i in 1:n_items) {
    connections <- all_connections(i, n_items)
    pairs <- connections_to_pair(connections)
    col_indices[, i] <- match(pairs, pair_names)
  }
  return(list(i = row_indices, j = c(col_indices), x = xes))
}

#' Get character representation of item connections
#' @param connections a data.frame returned by `all_connections`
connections_to_pair <- function(connections) {
  apply(connections, 1, function(x) paste0("x", paste0(x, collapse = "_")))
}

#' Find all connections of an element
#' @param i The point for which the connections are sought (must be in 1...n)
#' @param n The number of points
#' @return A `data.frame` with two columns representing start and end point.
#'   The second column always contains the "larger number" point.
all_connections <- function(i, n) {
  if (i > n | i <= 0)
    stop("Error in function `all_connections`: cannot find connections for element that is outside of legal range")
  connections <- expand.grid(i, setdiff(1:n, i))
  # put lower index in first column
  wrongorder <- connections[, 1] > connections[, 2]
  temp <- connections[, 2][wrongorder]
  connections[, 2][wrongorder] <- connections[, 1][wrongorder]
  connections[, 1][wrongorder] <- temp
  return(connections)
}
