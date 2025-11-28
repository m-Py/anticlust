
# Use blocking for anticlustering (anticluster within a block, which a level of a cagegorical variable, 
# or a stratum). When working on subsequent blocks, it "remembers" the assignment of prior blocks
# when optimizing an anticlustering objective.

# including constraints (must-link /cannot-link) would require re-indexing, I will at first not allow this
# because it is additional work for me.

blocked_anticlustering <- function(
    x, K, objective = "diversity", method = "exchange", 
    categories = NULL, cannot_link = NULL, blocks
) {
  stopifnot(method %in% c("exchange", "local-maximum"))

  validate_input(K, "K", len = 1) # here we can only generate equal-sized groups
  
  N <- nrow(x)

  if (!argument_exists(blocks)) stop("Can only use blocked anticlustering when argument 'blocks' is given.")  
  blocks <- merge_into_one_variable(blocks) # use "strata"
  
  n_blocks <- length(unique(blocks))
  blocksizes <- table(blocks)

  # Anticlustering with blocking: 
  condition_blocked <- rep(NA, N)
  # different initialization needed if there are cannot-link constraints
  if (argument_exists(cannot_link)) {
    condition_blocked <- optimal_cannot_link_reduced(
      N = N, K = K, 
      target_groups = table(initialize_clusters(N = N, K = K, NULL)),
      cannot_link = cannot_link
    )
    x <- convert_to_distances(x, squared = objective == "variance")
    x[cleanup_cannot_link_indices(cannot_link)] <- -(sum(x) + 1) # edit distances to maintain cannot-link restrictions
  }
  for (i in 1:n_blocks) {
    previous_groups <- condition_blocked # just for asserting at the end that no previous assignments were changed
    select <- blocks <= i
    if (is_distance_matrix(x)) {
      input <- x[select, select, drop = FALSE] 
    } else {
      input <- x[select, , drop = FALSE] 
    }
    target_groups <- table(initialize_clusters(N = sum(select), K = K, NULL))
    initial_groups <- add_unassigned_elements(target_groups, condition_blocked[select], N = N, K = K)
    # ensure that previous conditions are still as before:
    stopifnot(all(condition_blocked[select] == initial_groups, na.rm = TRUE))
    # Now do anticlustering with restriction that previously assigned conditions must not change.
    # Define who can be exchanged during anticlustering optimization (goes into the `categories` argument)
    exchange_partners <- get_blocking_exchange_parters(condition_blocked, select)
    condition_blocked[select] <- anticlustering(
      input, 
      K = initial_groups, 
      objective = objective,
      method = method,
      categories = cbind(exchange_partners, categories[select])
    )
    # ensure that previous conditions are still as before:
    stopifnot(all(condition_blocked[blocks < i] == previous_groups[blocks < i], na.rm = TRUE)) 
  }
  
  # verify that output has correct structure: 
  tab <- table(blocks, condition_blocked)
  stopifnot(all(abs(tab[ ,1] - tab[, 2]) <= 1)) # this test only works if equal group were requested, which are currently allowed only
  condition_blocked
}

# Create input for categories argument in anticlustering(), that ensures that previously blocked subjects
# are not changed (only used to compute the anticlustering objective)
get_blocking_exchange_parters <- function(condition_blocked, select) {
  exchange_partners <- rep(1, sum(select))
  is_already_assigned <- !is.na(condition_blocked[select])
  groups_already_assigned <- condition_blocked[select][is_already_assigned]
  exchange_partners[is_already_assigned] <- sample(sum(is_already_assigned), size = sum(is_already_assigned)) + 1
  exchange_partners
}