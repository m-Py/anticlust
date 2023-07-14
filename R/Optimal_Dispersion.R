
#' Maximize dispersion for K groups
#' 
#' @param x The data input. Can be one of two structures: (1) A
#'     feature matrix where rows correspond to elements and columns
#'     correspond to variables (a single numeric variable can be
#'     passed as a vector). (2) An N x N matrix dissimilarity matrix;
#'     can be an object of class \code{dist} (e.g., returned by
#'     \code{\link{dist}} or \code{\link{as.dist}}) or a \code{matrix}
#'     where the entries of the upper and lower triangular matrix
#'     represent pairwise dissimilarities.
#' @param K The number of groups or a vector describing the size of each group.
#' @param solver Optional argument; if passed, has to be either "glpk" or
#'   "symphony". See details.
#' @param max_dispersion_considered Optional argument used for early stopping. If the dispersion found
#'   is equal to or exceeds this value, a solution having the previous best dispersion 
#'   is returned.
#'
#' @export
#' 
#' @return A list with four elements:  
#'    \item{dispersion}{The optimal dispersion}
#'    \item{groups}{An assignment of elements to groups (vector)}
#'    \item{edges}{A matrix of 2 columns. Each row contains the indices of 
#'    elements that had to be investigated to find the dispersion (i.e., each pair
#'    of elements cannot be part of the same group in order to achieve maximum 
#'    dispersion).}
#'    \item{dispersions_considered}{All distances that were tested until the dispersion was found.}
#' 
#' @details
#'
#'   The dispersion is the minimum distance between two elements
#'   within the same group. This function implements an optimal method
#'   to maximize the dispersion. If the data input \code{x} is a feature
#'   matrix and not a dissimilarity matrix, the pairwise Euclidean
#'   distance is used. It uses the algorithm presented in Max
#'   Diekhoff's Bachelor thesis at the Computer Science Department at
#'   the Heinrich Heine University Düsseldorf.
#'
#'   To find out which items are not allowed to be grouped in the same
#'   cluster for maximum dispersion, the algorithm sequentially builds
#'   instances of a graph coloring problem, using an integer linear
#'   programming (ILP) representation (also see Fernandez et al.,
#'   2013).  It is possible to specify the ILP solver via the argument
#'   \code{solver}. This function either requires the R package
#'   \code{Rglpk} and the GNU linear programming kit
#'   (<http://www.gnu.org/software/glpk/>) or the R package
#'   \code{Rsymphony} and the COIN-OR SYMPHONY solver libraries
#'   (<https://github.com/coin-or/SYMPHONY>). If the argument
#'   \code{solver} is not specified, the function will try to find the
#'   GLPK or SYMPHONY solver by itself (it prioritizes using SYMPHONY if
#'   both are available). The GNU linear programming kit (\code{solver =
#'   "glpk"}) seems to be considerably slower for K >= 3 than the
#'   SYMPHPONY solver (\code{solver = "symphony"}).
#' 
#'   Optimally solving the maximum dispersion problem is NP-hard for K
#'   > 2 and therefore computationally infeasible for larger data
#'   sets. For K = 3 and K = 4, it seems that this approach scales up to several 100 elements, 
#'   or even > 1000 for K = 3 (at least when using the Symphony solver). 
#'   For larger data sets, use the heuristic approaches in \code{\link{anticlustering}} or
#'   \code{\link{bicriterion_anticlustering}}. However, note that for K = 2, 
#'   the optimal approach is usually much faster than the heuristics.
#'   
#'   In the output, the element \code{edges} defines which elements must be in separate 
#'   clusters in order to achieve maximum dispersion. All elements not listed here
#'   can be changed arbitrarily between clusters without reducing the dispersion.
#'   If the maximum possible dispersion corresponds to the minimum dispersion
#'   in the data set, the output elements \code{edges} and \code{groups} are set to
#'   \code{NULL} because all possible groupings have the same value of dispersion.
#'   In this case the output element \code{dispersions_considered} has length 1.
#'   
#'
#' @note If the SYMPHONY solver is used, an unfortunate
#' "message" is printed to the console when this function terminates: 
#' 
#' sym_get_col_solution(): No solution has been stored!
#' 
#' This message is no reason to worry and instead is a direct result
#' of the algorithm finding the optimal value for the dispersion.
#' Unfortunately, this message is generated in the C code underlying the 
#' SYMPHONY library (via the printing function \code{printf}), which cannot be
#' prevented in R.
#'
#' @author
#' 
#' Max Diekhoff
#' 
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @seealso \code{\link{dispersion_objective}} \code{\link{anticlustering}}
#' 
#' @references 
#' 
#' Diekhoff (2023). Maximizing dispersion for anticlustering. Retrieved from 
#' https://www.cs.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Algorithmische_Bioinformatik/Bachelor-_Masterarbeiten/2831963_ba_ifo_AbschlArbeit_klau_mapap102_madie120_20230203_1815.pdf  
#' 
#' Fernández, E., Kalcsics, J., & Nickel, S. (2013). The maximum dispersion 
#' problem. Omega, 41(4), 721–730. https://doi.org/10.1016/j.omega.2012.09.005 
#' 
#' 
#' @examples
#' 
#' N <- 30
#' M <- 5
#' K <- 3
#' data <- matrix(rnorm(N*M), ncol = M)
#' distances <- dist(data)
#' 
#' opt <- optimal_dispersion(distances, K = K)
#' opt
#' 
#' # Compare to bicriterion heuristic:
#' groups_heuristic <- anticlustering(
#'   distances, 
#'   K = K,
#'   method = "brusco", 
#'   objective = "dispersion", 
#'   repetitions = 100
#' )
#' c(
#'   OPT = dispersion_objective(distances, opt$groups),
#'   HEURISTIC = dispersion_objective(distances, groups_heuristic)
#' )
#' 
#' # Different group sizes are possible:
#' table(optimal_dispersion(distances, K = c(15, 10, 5))$groups)
#' 
#' # Induce cannot-link constraints by maximizing the dispersion:
#' solvable <- matrix(1, ncol = 6, nrow = 6)
#' solvable[2, 1] <- -1
#' solvable[3, 1] <- -1
#' solvable[4, 1] <- -1
#' solvable <- as.dist(solvable)
#' solvable
#' 
#' # An optimal solution has to put item 1 in a different group than 
#' # items 2, 3 and 4 -> this is possible for K = 2
#' optimal_dispersion(solvable, K = 2)$groups
#' 
#' # It no longer works when item 1 can also not be linked with item 5:
#' # (check out output!)
#' unsolvable <- as.matrix(solvable)
#' unsolvable[5, 1] <- -1
#' unsolvable <- as.dist(unsolvable)
#' unsolvable
#' optimal_dispersion(unsolvable, K = 2)
#' 
#'

optimal_dispersion <- function(x, K, solver = NULL, max_dispersion_considered = NULL) {
  
  if (argument_exists(solver)) {
    validate_input(solver, "solver", objmode = "character", len = 1,
                   input_set = c("glpk", "symphony"), not_na = TRUE, not_function = TRUE)
  } else {
    solver <- find_ilp_solver()
  }
  validate_data_matrix(x)
  validate_input(K, "K", objmode = "numeric", must_be_integer = TRUE, not_na = TRUE, not_function = TRUE)
  
  if (is.null(max_dispersion_considered)) {
    max_dispersion_considered <- Inf
  } else {
    validate_input(max_dispersion_considered, "max_dispersion_considered", 
                   objmode = "numeric", len = 1, not_na = TRUE, not_function = TRUE)
  }

  distances <- convert_to_distances(x)
  diag(distances) <- Inf
  N <- nrow(distances)
  if (!(length(K) == 1 || sum(K) == N)) {
    stop("Argument `K` is misspecified.")
  }
  
  # `target_groups` is primarily needed for unequal sized groups
  target_groups <- sort(table(initialize_clusters(N, K, NULL)), decreasing = TRUE)
  K <- length(target_groups)
  dispersion_found <- FALSE
  # Data frame to keep track of previous nearest neighbours (init as NULL)
  all_nns <- NULL
  # Placeholders to store data needed for retrieving the anticlusters
  last_solution <- NULL
  all_nns_last <- NULL
  all_nns_reordered_last <- NULL
  dispersions_considered <- NULL
  counter <- 1
  MINIMUM_DISTANCE <- min(distances)
  while (!dispersion_found) {
    dispersion <- min(distances)
    if (dispersion >= max_dispersion_considered) {
      break
    }
    ids_of_nearest_neighbours <- which(dispersion == distances, arr.ind = TRUE)
    all_nns <- rbind(all_nns, remove_redundant_edges(ids_of_nearest_neighbours))
    # Reorder edge labels so that they start from 1 to C, where C is the number
    # of relevant edges (Better for creating K-coloring ILP).
    all_nns_reordered <- reorder_edges(all_nns)
    # Construct graph from all previous edges (that had low distances)
    ilp <- k_coloring_ilp(all_nns_reordered, N, K, target_groups)
    solution <- solve_ilp_graph_colouring(ilp, solver)
    dispersion_found <- solution$status != 0 
    if (!dispersion_found){
      last_solution <- solution
      all_nns_last <- all_nns
      all_nns_reordered_last <- all_nns_reordered
      dispersions_considered <- c(dispersions_considered, dispersion)
    }
    counter <- counter + 1
    # Take out distances that have been investigated to proceed
    distances[ids_of_nearest_neighbours] <- Inf 
  }
  if (dispersion == MINIMUM_DISTANCE) { # no improvement for dispersion is possible; first test fails
    return(
      list(
        dispersion = MINIMUM_DISTANCE, 
        groups = NULL,
        edges = NULL, 
        dispersions_considered = MINIMUM_DISTANCE
      )
    )
  }
  # Calculate anticlusters from the last iteration with a K-coloring
  groups <- groups_from_k_coloring_mapping(
    result_value = last_solution$obj,
    result_x = last_solution$x,
    all_nns = all_nns_last,
    all_nns_reordered = all_nns_reordered_last,
    N = N,
    K = K, 
    target_groups = target_groups
  )
  return(
    list(
      dispersion = dispersion, 
      groups = groups,
      edges = unname(all_nns), # rownames can be quite ugly here
      dispersions_considered = c(dispersions_considered, dispersion)
    )
  )
}

k_coloring_ilp <- function(all_nns_reordered, N, K, target_groups) {
  # Initialize some constant variables
  nr_of_nodes <- max(all_nns_reordered)
  nr_of_edges <- nrow(all_nns_reordered)
  nr_of_x_variables <- nr_of_nodes * K
  equality_signs <- equality_identifiers()
  
  # Construct ILP constraint matrix
  constraints <- sparse_constraints_dispersion(all_nns_reordered, nr_of_nodes, nr_of_edges, K)
  colnames(constraints) <- constraint_names(nr_of_nodes, K)
  
  # Directions of the constraints:
  equalities <- c(rep(equality_signs$e, nr_of_nodes),
                  rep(equality_signs$l, (nr_of_edges +1) * K))
  
  # Right-hand-side of ILP
  rhs <- c(rep(1, nr_of_nodes), rep(0, nr_of_edges * K), target_groups)
  
  # Objective function of the ILP
  obj_function <- c(rep(1, K), rep(0, nr_of_x_variables))
  
  # Give names to all objects for inspection purposes
  names(obj_function) <- colnames(constraints)
  
  instance <- list()
  instance$constraints  <- constraints
  instance$equalities   <- equalities
  instance$rhs          <- rhs
  instance$obj_function <- obj_function
  
  return(instance)

}

# Construct a sparse matrix representing the ILP constraints
#
# @param all_nns_reordered A data frame where every row defines an edge in the graph
# @param nr_of_nodes The number of nodes in the graph
# @param K The maximum number of colors to be used
#
# @return A sparse matrix representing the left-hand side of the ILP (A in Ax ~ b)
#
sparse_constraints_dispersion <- function(all_nns_reordered, nr_of_nodes, nr_of_edges, K) {
  node <- node_constraints(nr_of_nodes, K)
  edge <- edge_constraints(all_nns_reordered, nr_of_nodes, nr_of_edges, K)
  group <- group_constraints_dispersion(nr_of_nodes, nr_of_edges, K)
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
group_constraints_dispersion <- function(nr_of_nodes, nr_of_edges, K)  {
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
groups_from_k_coloring_mapping <- function(result_value, result_x, all_nns, all_nns_reordered, N, K, target_groups) {
  
  # Retrieve assigned colors from the x_v,i ILP variables
  mapping <- rep(NA, max(all_nns_reordered))
  matrix <- matrix(result_x[-(1:K)], nrow=K)
  # sort labeling by frequency of occurrence
  matrix <- matrix[order(rowSums(matrix), decreasing = TRUE), ]
  for(l in 1:ncol(matrix)){
    mapping[l] <- which.max(matrix[,l]) # mapping is the groups, now retrieve original indices
  }
  
  df <- data.frame(
    mapping = mapping[c(all_nns_reordered)],
    tmp_id = c(all_nns_reordered),
    id = c(all_nns)
  )
  # create vector with groupings
  groups_new <- rep(NA, N)
  groups_new[df$id] <- df$mapping
  # now we have the original indices, assign remaining elements randomly to groups
  # how many are not yet assigned
  freq_not_assigned <- target_groups - table(groups_new)
  assigned <- table(groups_new)
  # Randomly fill other groups that are not yet full
  if (sum(assigned) < N) {
    groups_new[is.na(groups_new)] <- sample_(rep(1:K, freq_not_assigned))
  }
  groups_new
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
