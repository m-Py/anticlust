
#' Construct the ILP represenation of a cluster editing instance
#' 
#' (Can be used to solve anticlustering and clustering problems)
#'
#' @param distances An n x n matrix representing the
#'     distances between items
#' @param K The number of groups to be created
#' @param group_restriction If FALSE, there is no restriction on the
#'     group size, leading to a normal weighted cluster editing formulation
#'
#' @return A list representing the ILP formulation of the instance
#'
#' @noRd
#'

anticlustering_ilp <- function(distances, K, group_restriction = TRUE) {

  # Initialize some constant variables:
  equality_signs <- equality_identifiers()
  n        <- nrow(distances)
  group_size     <- n / K
  costs          <- vectorize_weights(distances)

  # Specify the number of triangular constraints:
  n_tris <- choose(n, 3) * 3

  # Construct ILP constraint matrix
  constraints <- sparse_constraints(n, costs$pair)
  colnames(constraints) <- costs$pair

  ## Directions of the constraints:
  equalities <- c(rep(equality_signs$l, n_tris),
                  rep(equality_signs$e, n))

  # Right-hand-side of ILP
  rhs <- c(rep(1, n_tris), rep(group_size - 1, n))

  # Objective function of the ILP
  obj_function <- costs$costs

  # Give names to all objects for inspection purposes
  names(obj_function) <- colnames(constraints)

  # For normal cluster editing, remove group constraints:
  if (group_restriction == FALSE) {
    rhs <- rhs[1:n_tris]
    equalities <- equalities[1:n_tris]
    constraints <- constraints[1:n_tris, ]
  }

  ## return instance
  instance              <- list()
  instance$n_groups     <- K
  instance$group_size   <- group_size
  instance$distances    <- distances
  instance$costs        <- costs
  instance$constraints  <- constraints
  instance$equalities   <- equalities
  instance$rhs          <- rhs
  instance$obj_function <- obj_function

  return(instance)
}

# Return identifiers for equality relationships
#
# @return A list of three elements containing strings representing
#     equality (e), lower (l), and greater (g) relationships
#
equality_identifiers <- function() {
  equal_sign <- "=="
  lower_sign <- "<="
  greater_sign <- ">="
  list(e = equal_sign, l = lower_sign, g = greater_sign)
}


# Convert matrix of distances into vector of distances
#
# @param distances A distance matrix of class matrix or dist
#
# @return A data.frame having the following columns: `costs` - the
#     distances in vectorized form; `i` the first index of the item
#     pair that is connected; `j` the second index of the item pair
#     that is connected; `pair` A string of form "xi_j" identifying the
#     item pair
vectorize_weights <- function(distances) {
  # Problem: I have matrix of costs but need vector for ILP.
  # Make vector of costs in data.frame (makes each cost identifiable)
  costs <- expand.grid(1:ncol(distances), 1:nrow(distances))
  colnames(costs) <- c("i", "j")
  costs$costs <- c(distances)
  # remove redundant or self distances:
  costs <- costs[costs$i < costs$j, ]
  costs$pair <- paste0("x", costs$i, "_", costs$j)
  rownames(costs) <- NULL
  return(costs)
}

# Construct a sparse matrix representing the ILP constraints
#
# @param n The number of elements
# @param pair_names A character vector of names representing the item
#     pairs (e.g.: "x1_2", ..., "x100_123"). Is read out from the data frame
#     have `costs` (i.e., costs$pair is passed as pair_names)
# @return A sparse matrix representing the left-hand side of the ILP (A
#     in Ax ~ b)
#
sparse_constraints <- function(n, pair_names) {
  tri <- triangular_constraints(n, pair_names)
  gr  <- group_constraints(n, pair_names)
  Matrix::sparseMatrix(c(tri$i, gr$i), c(tri$j, gr$j), x = c(tri$x, gr$x))
}

# Indices for sparse matrix representation of triangular constraints
triangular_constraints <- function(n, pair_names) {
  triangular_constraints <- choose(n, 3)
  coef_per_constraint <- 3
  # number of coefficients in constraint matrix:
  col_indices <- matrix(ncol = triangular_constraints, nrow = coef_per_constraint * 3)
  row_indices <- rep(1:(triangular_constraints*3), each = 3)
  xes <- rep(c(-1, 1, 1, 1, -1, 1, 1, 1, -1), triangular_constraints)
  # generate all item triplets
  triplets <- t(combn(1:n, 3))
  # get the respective column indices in constraint matrix
  ds <- c(t(data.frame(
    paste0("x", triplets[, 1], "_", triplets[, 2]),
    paste0("x", triplets[, 1], "_", triplets[, 3]),
    paste0("x", triplets[, 2], "_", triplets[, 3])
  )))
  col_indices <- matrix(match(ds, pair_names), nrow = 3)
  col_indices <- rbind(col_indices, col_indices, col_indices)
  list(i = row_indices, j = c(col_indices), x = xes)
}

# Indices for sparse matrix representation of group constraints
group_constraints <- function(n, pair_names) {
  coef_per_constraint <- (n - 1)
  group_coefficients <- coef_per_constraint * n
  row_indices <- rep((1:n) + (3 * choose(n, 3)), each = coef_per_constraint)
  col_indices <- matrix(ncol = n, nrow = coef_per_constraint)
  xes         <- rep(1, length = group_coefficients)
  for (i in 1:n) {
    connections <- all_connections(i, n)
    pairs <- paste0("x", connections[, 1], "_", connections[, 2])
    col_indices[, i] <- match(pairs, pair_names)
  }
  return(list(i = row_indices, j = c(col_indices), x = xes))
}

# Find all connections of an element
# @param i The point for which the connections are sought (must be in
#     1...n)
# @param n The number of points
# @return A `data.frame` with two columns representing start and end
#     point.  The second column always contains the "larger number"
#     point.
all_connections <- function(i, n) {
  connections <- expand.grid(i, setdiff(1:n, i))
  # put lower index in first column
  wrongorder <- connections[, 1] > connections[, 2]
  temp <- connections[, 2][wrongorder]
  connections[, 2][wrongorder] <- connections[, 1][wrongorder]
  connections[, 1][wrongorder] <- temp
  return(connections)
}
