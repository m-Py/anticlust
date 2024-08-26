
# Initialize must-link constraints by assigning all elements having the same ID to the same set
# all others remain free (i.e., as NA). This is a "randomized fit" algorithm for bin packing
get_init_assignments <- function(N, ID, target_groups) {
  # Initialize all as NA
  init <- rep(NA, N)
  K <- length(target_groups)
  cluster_sizes_real <- rep(0, K)
  multiple_IDs <- as.numeric(names(table(ID)[table(ID) > 1]))
  
  for (current_id in multiple_IDs) {
    random_order_clusters <- sample(K)
    for (k in random_order_clusters) {
      # only fill into cluster if it fits
      if ((cluster_sizes_real[k] + sum(ID == current_id)) > target_groups[k]) {
        if (k == random_order_clusters[K]) {
          stop("I could not fulfil the `must_link` restrictions, sorry!")
        }
        next
      }
      init[ID == current_id] <- k
      cluster_sizes_real[k] <- cluster_sizes_real[k] + sum(ID == current_id)
      break
    }
  }
  stopifnot(sum(!is.na(init))  == sum(ID %in% multiple_IDs))
  init
}

# Initialize the groupings for the reduced sample, for all elements:
init_must_link_groups <- function(N, IDs_initial, IDs_reduced, target_groups) {
  init <- get_init_assignments(N, IDs_initial, target_groups)
  init <- add_unassigned_elements(target_groups, init, N, length(target_groups))
  # only return one index per must-link group:
  init[sapply(IDs_reduced, FUN = "[", 1)]
}

# Adjust distances for must-link anticlustering (which uses a "reduced" data set where
# each must-link group is treated as a single unit)
adjusted_distances_must_link <- function(distances, must_link) {
  distances <- as.matrix(distances)
  N <- nrow(distances)
  stopifnot(N == length(must_link))
  N_reduced <- length(unique(must_link))
  new_distances <- matrix(NA, ncol = N_reduced, nrow = N_reduced)
  list_must_link_indices <- tapply(1:N, must_link, c)
  for (i in 1:N_reduced) {
    for (j in 1:N_reduced) {
      new_distances[i, j] <- sum(distances[list_must_link_indices[[i]], list_must_link_indices[[j]]])
    }
  }
  list(distances = new_distances, IDs = list_must_link_indices)
}

#' This is the function that is called from anticlustering()
must_link_anticlustering <- function(x, K, must_link, method = "exchange", objective = "diversity", repetitions = NULL) {
  validate_input(objective, "objective", input_set = c("diversity"), not_na = TRUE, len = 1) 
  validate_input(method, "method", input_set = c("local-maximum", "exchange"), not_na = TRUE, len = 1) 
  stopifnot(is_distance_matrix(x))
  x <- to_matrix(x)
  N <- nrow(x)
  
  must_link <- replace_na_by_index(must_link)
  
  dt <- adjusted_distances_must_link(distances, must_link)
  
  target_groups <- table(initialize_clusters(N, K, NULL))
  
  # possibly use multiple initializing partitions:
  if (argument_exists(repetitions)) {
    init_partitions <- t(replicate(
      n = repetitions,
      init_must_link_groups(N, IDs_initial = must_link, IDs_reduced = dt$IDs, target_groups = target_groups)
    ))
    init_partitions <- init_partitions - 1 # yeah, this is not good code (where should this go?)
    init <- init_partitions[1, ]
  } else {
    init_partitions <- NULL
    init <- init_must_link_groups(N, IDs_initial = must_link, IDs_reduced = dt$IDs, target_groups = target_groups)
  }
  
  reduced_clusters <- c_anticlustering(
    dt$distances,
    K = init,
    categories = lengths(dt$IDs), # restrict exchanges to node with same number of elements
    objective = "diversity",
    local_maximum = ifelse(method == "local-maximum", TRUE, FALSE),
    init_partitions = init_partitions
  )
  
  full_clusters <- rep(NA, N)
  for (i in seq_along(dt$IDs)) {
    full_clusters[dt$IDs[[i]]] <- reduced_clusters[i]
  }
  full_clusters
}
