# Optimal algorithm for one dimensional bin packing

optimal_binpacking_ <- function(capacities, weights, solver = NULL, time_limit = NULL) {
  
  n_batches <- length(capacities)
  n <- length(weights)
  
  # We get names for the decision variables to keep track of everything
  variables <- get_col_names(n_batches, n)
  
  # Decision variables: Is item i assigned to batch j (so there are 
  # N x C constraints for N items and C batches)
  
  # Two types of constraints:
  # a. Capacity of batches must not be exceeded (1 constraint per batch)
  constraints1 <- matrix(0, ncol = length(variables), nrow = n_batches)
  colnames(constraints1) <- variables
  for (i in 1:nrow(constraints1)) {
    constraints1[i, grepl(paste0("b", i, "_"), variables)] <- weights
  }
  
  # b. Each item is filled into exactly one bin (1 constraint per item)
  constraints2 <- matrix(0, ncol = length(variables), nrow = n)
  colnames(constraints2) <- variables
  for (i in 1:nrow(constraints2)) {
    constraints2[i, grepl(paste0("item_", i, "_"), variables)] <- 1
  }
  
  dirs <- c(rep("<=", nrow(constraints1)), rep("==", nrow(constraints2)))
  rhs <- c(capacities, rep(1, nrow(constraints2)))
  
  ilp <- list()
  ilp$constraints  <- rbind(constraints1, constraints2)
  ilp$equalities   <- dirs
  ilp$rhs          <- rhs
  ilp$obj_function <- rep(1, ncol(constraints1))
  
  ilp_solution <- solve_ilp(ilp, objective = "min", solver = solver, time_limit = time_limit)
  if (ilp_solution$status != 0) {
    stop("The constraints cannot be fulfilled (really).")
  }
  target_bins <- variables[ilp_solution$x == 1]
  # we have to extract the actual target bin from the strings, like this:
  as.integer(gsub("b", "", sapply(strsplit(target_bins, "_"), "[", 1)))
}

get_col_names <- function(n_batches, n) {
  combs <- expand.grid(1:n_batches, 1:n)
  paste0("b", combs[,1], "_item_", combs[, 2], "_")
}
