
#' Maximize dispersion for K groups
#' 
#' @param distances Input data, dissimilarity matrix
#' @param K The number of anticlusters to be created
#' @param solver Either "glpk" (default) or "symphony"
#'
#' @export
#' 
#' @return A list with three elements:  
#'    \code{dispersion}: The optimal dispersion; 
#'    \code{groups}: A possible assignment to groups (vector);
#'    \code{edges}: A matrix of 2 columns. Each row contains the indices of 
#'    elements that had to be investigated to find the dispersion. (i.e., each pair
#'    of elements cannot be part of the same group in order to achieve maximum 
#'    dispersion)
#' 
#' @details Finds the optimal dispersion as using the algorithm presented in 
#'   Max Diekhoff's Bachelor thesis. 
#' 
#' @author
#' 
#' Max Diekhoff
#' 
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @references 
#' 
#' Diekhoff (2023). Maximizing dispersion for anticlustering. Retrieved from https://www.cs.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Algorithmische_Bioinformatik/Bachelor-_Masterarbeiten/2831963_ba_ifo_AbschlArbeit_klau_mapap102_madie120_20230203_1815.pdf  
#' 
#' Fernández, E., Kalcsics, J., & Nickel, S. (2013). The maximum dispersion problem. Omega, 41(4), 721–730. https://doi.org/10.1016/j.omega.2012.09.005 
#'

optimal_dispersion <- function(distances, K, solver = "glpk") {
  
  validate_input(solver, "solver", objmode = "character", len = 1,
                 input_set = c("glpk", "symphony"), not_na = TRUE, not_function = TRUE)
  
  distances <- as.matrix(distances)
  diag(distances) <- Inf
  N <- nrow(distances)
  dispersion_found <- FALSE
  # Data frame to keep track of previous nearest neighbours (init as NULL)
  all_nns <- NULL
  # Placeholders to store data needed for retrieving the anticlusters
  last_solution <- NULL
  all_nns_last <- NULL
  all_nns_reordered_last <- NULL
  counter <- 1
  # TODO: Include code on what to do when the first check for bipartiteness already fails
  while (!dispersion_found) {
    dispersion <- min(distances)
    ids_of_nearest_neighbours <- which(dispersion == distances, arr.ind = TRUE)
    all_nns <- rbind(all_nns, remove_redundant_edges(ids_of_nearest_neighbours))
    # Reorder edge labels so that they start from 1 to C, where C is the number
    # of relevant edges (Better for creating K-coloring ILP).
    all_nns_reordered <- reorder_edges(all_nns)
    # Construct graph from all previous edges (that had low distances)
    #gr <- graph_from_edgelist(all_nns_reordered, directed = FALSE)
    ilp <- k_coloring_ilp(all_nns_reordered, N, K)
    solution <- solve_ilp_graph_colouring(ilp, solver)
    # Dispersion is found when graph has no k-coloring, which corresponds to objective value beeing 0
    dispersion_found <- solution$status != 0
    if(!dispersion_found){
      last_solution <- solution
      all_nns_last <- all_nns
      all_nns_reordered_last <- all_nns_reordered
    }
    counter <- counter +1
    # Take out distances that have been investigated to proceed
    distances[ids_of_nearest_neighbours] <- Inf 
  }
  # Calculate anticlusters from the last iteration with a K-coloring
  groups <- groups_from_k_coloring_mapping(
    last_solution$obj,
    last_solution$x,
    all_nns_last,
    all_nns_reordered_last,
    N,
    K
  )
  return(
    list(
      dispersion = dispersion, 
      groups = groups$groups,
      edges = all_nns
    )
  )
}

k_coloring_ilp <- function(all_nns_reordered, N, K){
  # Initialize some constant variables
  nr_of_nodes <- max(all_nns_reordered)
  nr_of_edges <- nrow(all_nns_reordered)
  nr_of_x_variables <- nr_of_nodes * K
  max_group_size <- ceiling(N / K)
  equality_signs <- equality_identifiers()
  
  # Construct ILP constraint matrix
  constraints <- sparse_constraints_dispersion(all_nns_reordered, nr_of_nodes, nr_of_edges, max_group_size, K)
  colnames(constraints) <- constraint_names(nr_of_nodes, K)
  
  # Directions of the constraints:
  equalities <- c(rep(equality_signs$e, nr_of_nodes),
                  rep(equality_signs$l, (nr_of_edges +1) * K))
  
  # Right-hand-side of ILP
  rhs <- c(rep(1, nr_of_nodes), rep(0, nr_of_edges * K), rep(max_group_size, K))
  
  # Objective function of the ILP
  obj_function <- c(rep(1, K), rep(0, nr_of_x_variables))
  
  # Give names to all objects for inspection purposes
  names(obj_function) <- colnames(constraints)
  
  instance <- list()
  instance$n_groups     <- K
  instance$group_size   <- max_group_size
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

# Construct a sparse matrix representing the ILP constraints
#
# @param all_nns_reordered A data frame where every row defines an edge in the graph
# @param nr_of_nodes The number of nodes in the graph
# @param max_group_size The maximum group size which corresponds the the maximum nr of nodes that are allowed to be assigned the same color
# @param K The maximum number of colors to be used
#
# @return A sparse matrix representing the left-hand side of the ILP (A in Ax ~ b)
#
sparse_constraints_dispersion <- function(all_nns_reordered, nr_of_nodes, nr_of_edges, max_group_size, K) {
  node <- node_constraints(nr_of_nodes, K)
  edge <- edge_constraints(all_nns_reordered, nr_of_nodes, nr_of_edges, K)
  group <- group_constraints_dispersion(nr_of_nodes, nr_of_edges, max_group_size, K)
  Matrix::sparseMatrix(c(node$i, edge$i, group$i), c(node$j, edge$j, group$j), x = c(node$x, edge$x, group$x))
}

# Indices for sparse matrix representation of node constraints
node_constraints <- function(nr_of_nodes, K) {
  col_indices <- rep(1:nr_of_nodes, each=K)
  row_start <- K+1
  row_end <- (nr_of_nodes+1)*K
  row_indices <- row_start:row_end
  xes <- rep(1, nr_of_nodes * K)
  list(i = col_indices, j = row_indices, x = xes)
}

# Indices for sparse matrix representation of edge constraints
edge_constraints <- function(all_nns, nr_of_nodes, nr_of_edges, K) {
  col_start <- nr_of_nodes+1
  col_end <- nr_of_nodes+(nr_of_edges * K)
  col_indices <- rep(col_start:col_end, each=3)
  row_indices <- c()
  for(i in 1:nr_of_edges){
    u <- all_nns[i,1]
    v <- all_nns[i,2]
    current_rows <- c(rep(c(0, u * K, v * K), K))
    row_indices <- c(row_indices, current_rows + rep(1:K, each=3))
  }
  xes <- rep(c(-1,1,1), nr_of_edges * K)
  list(i = col_indices, j = row_indices, x = xes)
}

# Indices for sparse matrix representation of group size constraints
group_constraints_dispersion <- function(nr_of_nodes, nr_of_edges, max_group_size, K)  {
  col_start <- nr_of_nodes+(nr_of_edges * K)+1
  col_end <- nr_of_nodes+(nr_of_edges * K)+K
  col_indices <- rep(col_start:col_end, each=nr_of_nodes)
  row_indices <- (rep(1:nr_of_nodes) * K) + rep(1:K, each=nr_of_nodes) 
  xes <- rep(1, nr_of_nodes * K)
  list(i = col_indices, j = row_indices, x = xes)
}


# Return variable names to have meaningful column names for the objective function and constraint matrix
#
# @param nr_of_nodes The number of nodes in the graph
# @param K The maximum number of colors to be used
#
# @return A vector with w_l (l=1,...,K) and x_i_j (i=1,...,nr_of_nodes, j=1,...,K)
#
constraint_names <- function(nr_of_nodes, K) {
  w_l <- rep("w_1", 1)
  for(l in 2:K){
    w_l <- c(w_l, paste0("w_", l))
  }
  
  x_j_i <- expand.grid(1:K, 1:nr_of_nodes)
  colnames(x_j_i) <- c("j", "i")
  x_j_i$variables <- paste0("x_",x_j_i$i,"_",x_j_i$j)
  
  return(c(w_l, x_j_i$variables))
}

# allow for different solvers (Symphony, GLPK)

solve_ilp_graph_colouring <- function(ilp, solver) {
  
  # solver_function = Rglpk::Rglpk_solve_LP OR Rsymphony::Rsymphony_solve_LP
  # name_opt (refers to the output of the function) = "objval" (GLPK) OR "optimum" (Symphony)
  # rest of the input is the same between the solver functions, which is nice
  solver_function <- ifelse(solver == "symphony", Rsymphony::Rsymphony_solve_LP, Rglpk::Rglpk_solve_LP)
  name_opt <- ifelse(solver == "symphony", "objval", "optimum")
  
  ilp_solution <- solver_function(
    obj = ilp$obj_function,
    mat = ilp$constraints,
    dir = ilp$equalities,
    rhs = ilp$rhs,
    types = "B",
    max = FALSE
  )

  # return the optimal value and the variable assignment
  ret_list <- list() 
  ret_list$x <- ilp_solution$solution
  ret_list$obj <- ilp_solution[[name_opt]]
  ret_list$status <- ilp_solution$status
  ## name the decision variables
  names(ret_list$x) <- colnames(ilp$constraints)
  ret_list
}


# Restore a grouping from the solved ilp
groups_from_k_coloring_mapping <- function(result_value, result_x, all_nns, all_nns_reordered, N, K) {
  # Retrieve assigned colors from the x_v,i ILP variables
  mapping <- rep(NA, max(all_nns_reordered))
  matrix <- matrix(result_x[-(1:K)], nrow=K)
  for(l in 1:ncol(matrix)){
    mapping[l] <- which.max(matrix[,l])
  }
  
  df <- data.frame(
    mapping = mapping[c(all_nns_reordered)],
    tmp_id = c(all_nns_reordered),
    id = c(all_nns)
  )
  # create vector with groupings
  groups_new <- rep(NA, N)
  groups_new[df$id] <- df$mapping
  # how many are not yet assigned
  freq_not_assigned <- rep(N / K, K)
  assigned <- table(groups_new)
  for(i in 1:result_value){
    freq_not_assigned[i] <- freq_not_assigned[i] - assigned[i]
  }
  groups_new[is.na(groups_new)] <- sample(rep(1:K, freq_not_assigned))
  list(groups = groups_new, edges = df)
}

# Function that takes an edge list with any integer indices and remaps
# these indices to 1, ..., C (where C is the number of nodes) while 
# preserving the original order of indices. This is important for igraph
reorder_edges <- function(edgelist) {
  dims <- dim(edgelist)
  edgelist <- as.numeric(as.factor(edgelist))
  dim(edgelist) <- dims
  edgelist
}

remove_redundant_edges <- function(df) {
  df <- t(apply(df, 1, sort))
  df[!duplicated(df), ]
}
