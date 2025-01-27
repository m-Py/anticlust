

## The functions `get_init_assignments*()` get the assignments for the items
# that are restricted via cannot-link constraints. The other items still need 
# to be assigned (via `add_unassigned_elements()`, which is done in 
# init_must_link_groups()). 
# - get_init_assignments_optimal() implements an optimal
# bin packing algorithm. It is just a decision version because we are not interested
# in minimizing the number of bins (so the objective function is constant for 
# all possible assignments). 
# - get_init_assignments_heuristic() implements a randomized fit heuristic,
# where each must-link cluster is assigned to a random group where it fits.
# Only if this heuristic fails, the optimal algorithm will be called.

get_init_assignments <- function(N, ID, target_groups, method = "heuristic") {
  if (method == "optimal") {
    get_init_assignments_optimal(N, ID, target_groups)
  } 
  get_init_assignments_heuristic(N, ID, target_groups)
}

get_init_assignments_optimal <- function(N, ID, target_groups) {
  weights <- table(ID)[table(ID) > 1]
  multiple_IDs <- as.numeric(names(weights))
  opt_assignment <- optimal_binpacking_(target_groups, weights)
  # Optimal assignment is done for the "reduced" data set, we have to assign cluster
  # label to all elements that are part of must-link group:
  full_clusters <- rep(NA, N)
  for (i in seq_along(multiple_IDs)) {
    full_clusters[ID == multiple_IDs[i]] <- opt_assignment[i]
  }
  full_clusters
}

# Initialize must-link constraints by assigning all elements having the same ID to the same set
# all others remain free (i.e., as NA). This is a "randomized fit" algorithm for bin packing
get_init_assignments_heuristic <- function(N, ID, target_groups) {
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
  init <- tryCatch(
    get_init_assignments(N, IDs_initial, target_groups, method = "heuristic"),
    error = function(e) e
  )
  if ("simpleError" %in% class(init)) {
    init <- tryCatch(
      get_init_assignments(N, IDs_initial, target_groups, method = "optimal"),
      error = function(e) e
    )
  }
  if ("simpleError" %in% class(init)) {
    stop("The must-link constraints cannot be fulfilled! I really tried.")
  }
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
  list_must_link_indices <- get_must_link_indices(must_link)
  for (i in 1:N_reduced) {
    for (j in 1:N_reduced) {
      new_distances[i, j] <- sum(distances[list_must_link_indices[[i]], list_must_link_indices[[j]]])
    }
  }
  list(distances = new_distances, IDs = list_must_link_indices)
}

# This is the function that is called from anticlustering()
must_link_anticlustering <- function(x, K, must_link, method = "exchange", objective = "diversity", repetitions = NULL) {
  
  x <- to_matrix(x)
  N <- nrow(x)
  stopifnot(is_distance_matrix(x)) # caller must ensure that the distance matrix is here
  
  must_link <- to_numeric(must_link)
  must_link <- replace_na_by_index(must_link)
  
  dt <- adjusted_distances_must_link(x, must_link)
  
  target_groups <- table(initialize_clusters(N, K, NULL))
  
  # possibly use multiple initializing partitions:
  if (argument_exists(repetitions)) {
    init_partitions <- t(replicate(
      n = ifelse(method == "2PML", max(round(repetitions/2), 1), repetitions), # half repetitions for first phase
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
    local_maximum = ifelse(method == "exchange", FALSE, TRUE),
    init_partitions = init_partitions
  )
  
  full_clusters <- merged_cluster_to_original_cluster(reduced_clusters, must_link)
  
  ## This is the second phase of the 2PML algorithm:
  if (method == "2PML") {
    iterations <- ifelse(is.null(repetitions), 1, max(round(repetitions/2), 1)) # half repetitions for second phase
    BEST <- diversity_objective_(full_clusters, x)
    for (i in 1:iterations) {
      full_clusters_new <- improvement_search_2pml(x, full_clusters, must_link)
      # full clusters new is a perturbed partition, not locally optimal -> restore local optimality!
      reduced_clusters_new <- original_cluster_to_merged_cluster(full_clusters_new, must_link)
      # restore local optimality
      reduced_clusters_new <- c_anticlustering(
        dt$distances,
        K = reduced_clusters_new,
        categories = lengths(dt$IDs), # restrict exchanges to node with same number of elements
        objective = "diversity",
        local_maximum = TRUE
      )
      full_clusters <- merged_cluster_to_original_cluster(reduced_clusters_new, must_link)
    }
  }
  full_clusters
}

# Given an assignment of original elements to clusters, generate the
# corresponding clusters for the merged elements. 
original_cluster_to_merged_cluster <- function(clusters, must_link) {
  df <- data.frame(clusters, must_link)
  new_df <- df[!duplicated(df$must_link), ]
  # order by must_link grouping
  new_df[order(new_df$must_link), ]$clusters
}

# Given an assignment of "merged" elements to clusters, re-establish the
# corresponding clusters for the original elements. 
merged_cluster_to_original_cluster <- function(merged_clusters, must_link) {
  df <- data.frame(must_link, order = 1:length(must_link))
  new_order <- order(must_link)
  df <- df[new_order, ]
  df$clusters <- rep(merged_clusters, table(must_link))
  df[order(df$order), "clusters"]
}

## Do an improvement phase to overcome bad local optima
# For each must-link clique
improvement_search_2pml <- function(x, full_clusters, must_link) {
  N <- length(full_clusters)
  # get current objective for optimization
  OBJ_BY_CLUSTER <- diversity_objective_by_group(full_clusters, x)
  OBJ <- sum(OBJ_BY_CLUSTER)
  ## after this: set diagonal to zero, for more convenient updates of objective
  diag(x) <- 0

  # book keeping of indices (separately for must-linked elements and "singletons" that are not must-linked)
  list_must_link_indices <- get_must_link_indices(must_link) 
  cliques <- get_cliques(list_must_link_indices) 
  singletons <- get_singletons(list_must_link_indices)
  singleton_clusters <- full_clusters[singletons]
  
  # do swapping phase; iterate over cliques
  # for each clique, generate exchange partners twice
  for (i in seq_along(cliques)) {
    n_clique <- length(cliques[[i]])
    # randomly determine if the clique is exchanged with singletons only or with other cliques (that is, clique(s) + possibly singleton(s) in combination) 
    singleton_change <- sample(c(TRUE, FALSE), size = 1)
    if (singleton_change || n_clique <= 2) { # for cliques of size 2, we can only exchange with two singletons
      exchange_cluster <- get_exchange_partners_singletons( # generate exchange partner from singletons / not part of clique
        singletons, 
        singleton_clusters, 
        n_clique
      )
    } else { # only do the clique swapping for clique larger than 2 samples
      exchange_cluster <- get_exchange_partners_clique(cliques, i, full_clusters, must_link) # generate exchange partner from other cliques
    }
    if (!is.null(exchange_cluster$cluster_id)) {
      ## Do the swap
      tmp_clusters <- full_clusters
      tmp <- tmp_clusters[cliques[[i]]][1]
      tmp_clusters[cliques[[i]]] <- exchange_cluster$cluster_id
      tmp_clusters[exchange_cluster$sample_ids] <- tmp
      # perform swap if it improves objective
      ## Use local updating here instead of recomputing entirely, improves speed quite a bit
      OBJ_BY_CLUSTER_NEW <- update_diversity_must_link(
        x, OBJ_BY_CLUSTER, 
        full_clusters, 
        tmp_clusters, 
        indices_1 = cliques[[i]], 
        indices_2 = exchange_cluster$sample_ids, 
        cluster_1 = tmp, 
        cluster_2 = exchange_cluster$cluster_id
      )
      OBJ_NEW <- sum(OBJ_BY_CLUSTER_NEW)
      if (OBJ_NEW > OBJ) {
        OBJ <- OBJ_NEW
        OBJ_BY_CLUSTER <- OBJ_BY_CLUSTER_NEW
        full_clusters <- tmp_clusters
        singleton_clusters <- full_clusters[singletons]
      }
    }
  }
  full_clusters
}

update_diversity_must_link <- function(x, OBJ_BY_CLUSTER, clusters_old, clusters_new, indices_1, indices_2, cluster_1, cluster_2) {
  sum_distances_element_1_cluster_1 <- sum(x[indices_1, , drop = FALSE][, clusters_old == cluster_1, drop = FALSE])
  sum_distances_element_1_cluster_2 <- sum(x[indices_1, , drop = FALSE][, clusters_new == cluster_2, drop = FALSE])
  sum_distances_element_2_cluster_2 <- sum(x[indices_2, , drop = FALSE][, clusters_old == cluster_2, drop = FALSE])
  sum_distances_element_2_cluster_1 <- sum(x[indices_2, , drop = FALSE][, clusters_new == cluster_1, drop = FALSE])
  OBJ_BY_CLUSTER[cluster_1] <- OBJ_BY_CLUSTER[cluster_1] + sum_distances_element_2_cluster_1 - sum_distances_element_1_cluster_1
  OBJ_BY_CLUSTER[cluster_2] <- OBJ_BY_CLUSTER[cluster_2] + sum_distances_element_1_cluster_2 - sum_distances_element_2_cluster_2
  OBJ_BY_CLUSTER
}

# get list of indices of must-link cliques (and singletons)
get_must_link_indices <- function(must_link) {
  N <- length(must_link)
  tapply(1:N, must_link, c)
}

# using the list of must-link indices as input, get singleton indices (as vector)
get_singletons <- function(list_must_link_indices) {
  unname(unlist(list_must_link_indices[lengths(list_must_link_indices) == 1]))
}

#using the list  of must-link indices as input, get indices of cliques  (as list)
get_cliques <- function(list_must_link_indices) {
  list_must_link_indices[lengths(list_must_link_indices) > 1]
}

# determine a random cluster that fits a clique (and sample IDs for the exchange)
get_exchange_partners_singletons <- function(singletons, singleton_clusters, n_clique) {
  tab <- table(singleton_clusters)
  fitting_clusters <- as.numeric(names(tab[tab >= n_clique]))
  # If no swap attempt is possible, return list will NULL
  if (length(fitting_clusters) == 0) {
    return(list(
        cluster_id = NULL, 
        sample_ids = NULL
    ))
  }
  cluster_id <- sample_(fitting_clusters, size = 1)
  list(
    cluster_id = cluster_id, 
    sample_ids = sample(singletons[singleton_clusters == cluster_id], size = n_clique)
  )
}


## This function is quite horrible right now and I hope to make it better
# It generates multiple exchange partners for a must-link clique, using
# a combination of cliques / singletons. It randomly selects one combination that
# fits. First it has to generate all possible combinations that generates 
# the size of a must-link group (e.g., 6 = 1+2+3; 2+4; 5+1). This is quite hard
# in general I think.
get_exchange_partners_clique <- function(cliques, index, full_clusters, must_link) {
  
  n_clique <- length(cliques[[index]])
  stopifnot(n_clique > 2) # this must be here

  nl <- valid_sums_clique(n_clique) # returns a list of combinations that sum to n_clique
  
  # now check if any of these combinations exists in a cluster
  tab <- table(full_clusters, must_link)
  # iterate through combinations in random order and use the first exchange partners that work
  random_order <- sample(length(nl))
  for (i in random_order) {
    # determine where all samples are available
    cluster_ids_that_fit <- which(colSums(!apply(tab, 1, function(x) nl[[i]] %in% x)) == 0)
    # found fit!
    if (length(cluster_ids_that_fit) > 0) {
      one_cluster_id_that_fit <- sample_(cluster_ids_that_fit, size = 1)
      selected_cluster <- full_clusters == one_cluster_id_that_fit
      # select IDs of elements that serve as exchange ṕartners, via size of cliques
      tab2 <- table(full_clusters[selected_cluster], must_link[selected_cluster])
      must_link_ids <- as.numeric(dimnames(tab2)[[2]][random_match(nl[[i]], tab2)])
      # from IDs from selected must-link groups, get sample IDs
      sample_ids <- which(must_link %in% must_link_ids)
      if (any(is.na(must_link_ids)) || length(sample_ids) != n_clique) {
        next
      }
      return(list(
        cluster_id = one_cluster_id_that_fit, 
        sample_ids = sample_ids
      ))
    }

  }
  return(list(
    cluster_id = NULL, 
    sample_ids = NULL
  ))
}

## I generated these sums separately because recomputation is too expensive
valid_sums_clique <- function(n_clique) { # return a list of combinations that sum to n_clique
  if (n_clique <= 10) {
    return(
        list("3" = list(1:2), "4" = list(c(1L, 3L), c(2L, 2L), c(1L, 
1L, 2L)), "5" = list(c(1L, 4L), 2:3, c(1L, 1L, 3L), c(1L, 2L, 
2L), c(1L, 1L, 1L, 2L)), "6" = list(c(1L, 5L), c(2L, 4L), c(3L, 
3L), c(1L, 1L, 4L), 1:3, c(2L, 2L, 2L), c(1L, 1L, 1L, 3L), c(1L, 
1L, 2L, 2L), c(1L, 1L, 1L, 1L, 2L)), "7" = list(c(1L, 6L), c(2L, 
5L), 3:4, c(1L, 1L, 5L), c(1L, 2L, 4L), c(1L, 3L, 3L), c(2L, 
2L, 3L), c(1L, 1L, 1L, 4L), c(1L, 1L, 2L, 3L), c(1L, 2L, 2L, 
2L), c(1L, 1L, 1L, 1L, 3L), c(1L, 1L, 1L, 2L, 2L), c(1L, 1L, 
1L, 1L, 1L, 2L)), "8" = list(c(1L, 7L), c(2L, 6L), c(3L, 5L), 
    c(4L, 4L), c(1L, 1L, 6L), c(1L, 2L, 5L), c(1L, 3L, 4L), c(2L, 
    2L, 4L), c(2L, 3L, 3L), c(1L, 1L, 1L, 5L), c(1L, 1L, 2L, 
    4L), c(1L, 1L, 3L, 3L), c(1L, 2L, 2L, 3L), c(2L, 2L, 2L, 
    2L), c(1L, 1L, 1L, 1L, 4L), c(1L, 1L, 1L, 2L, 3L), c(1L, 
    1L, 2L, 2L, 2L), c(1L, 1L, 1L, 1L, 1L, 3L), c(1L, 1L, 1L, 
    1L, 2L, 2L), c(1L, 1L, 1L, 1L, 1L, 1L, 2L)), "9" = list(c(1L, 
8L), c(2L, 7L), c(3L, 6L), 4:5, c(1L, 1L, 7L), c(1L, 2L, 6L), 
    c(1L, 3L, 5L), c(1L, 4L, 4L), c(2L, 2L, 5L), 2:4, c(3L, 3L, 
    3L), c(1L, 1L, 1L, 6L), c(1L, 1L, 2L, 5L), c(1L, 1L, 3L, 
    4L), c(1L, 2L, 2L, 4L), c(1L, 2L, 3L, 3L), c(2L, 2L, 2L, 
    3L), c(1L, 1L, 1L, 1L, 5L), c(1L, 1L, 1L, 2L, 4L), c(1L, 
    1L, 1L, 3L, 3L), c(1L, 1L, 2L, 2L, 3L), c(1L, 2L, 2L, 2L, 
    2L), c(1L, 1L, 1L, 1L, 1L, 4L), c(1L, 1L, 1L, 1L, 2L, 3L), 
    c(1L, 1L, 1L, 2L, 2L, 2L), c(1L, 1L, 1L, 1L, 1L, 1L, 3L), 
    c(1L, 1L, 1L, 1L, 1L, 2L, 2L), c(1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 2L)), "10" = list(c(1L, 9L), c(2L, 8L), c(3L, 7L), c(4L, 
6L), c(5L, 5L), c(1L, 1L, 8L), c(1L, 2L, 7L), c(1L, 3L, 6L), 
    c(1L, 4L, 5L), c(2L, 2L, 6L), c(2L, 3L, 5L), c(2L, 4L, 4L
    ), c(3L, 3L, 4L), c(1L, 1L, 1L, 7L), c(1L, 1L, 2L, 6L), c(1L, 
    1L, 3L, 5L), c(1L, 1L, 4L, 4L), c(1L, 2L, 2L, 5L), 1:4, c(1L, 
    3L, 3L, 3L), c(2L, 2L, 2L, 4L), c(2L, 2L, 3L, 3L), c(1L, 
    1L, 1L, 1L, 6L), c(1L, 1L, 1L, 2L, 5L), c(1L, 1L, 1L, 3L, 
    4L), c(1L, 1L, 2L, 2L, 4L), c(1L, 1L, 2L, 3L, 3L), c(1L, 
    2L, 2L, 2L, 3L), c(2L, 2L, 2L, 2L, 2L), c(1L, 1L, 1L, 1L, 
    1L, 5L), c(1L, 1L, 1L, 1L, 2L, 4L), c(1L, 1L, 1L, 1L, 3L, 
    3L), c(1L, 1L, 1L, 2L, 2L, 3L), c(1L, 1L, 2L, 2L, 2L, 2L), 
    c(1L, 1L, 1L, 1L, 1L, 1L, 4L), c(1L, 1L, 1L, 1L, 1L, 2L, 
    3L), c(1L, 1L, 1L, 1L, 2L, 2L, 2L), c(1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 3L), c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L), c(1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 2L)))[[as.character(n_clique)]]
    )
  } else {
    # heuristic: return equal-sized cliques, but add 1s to compensate if n_clique is not divisible by X
    # (e.g., for n_clique = 12, we have 6x2, 4x3, 2x5+2x1, 2x6)
    max_number <- floor(n_clique/2)
    candidates <- get_set_for_subset_problem(n_clique)
    candidates <- candidates[candidates <= max_number] 
    subsets_init <- lapply(1:max_number, function(x) candidates[candidates == x])[-1] # first entry is the 1s, which we don't use here
    diffs <- lapply(lapply(subsets_init, sum), "-", n_clique) # get difference of sum to each 
    lapply(subsets_init, function(x) c(x, rep(1, n_clique - sum(x))))
  }
}

get_set_for_subset_problem <- function(n_clique) {
  init <- rep(1:n_clique, n_clique:1)
  init_list <- mapply(FUN = rep, 1:n_clique, n_clique:1)
  cumsums <- lapply(1:n_clique, function(x) cumsum(init[init == x]))
  valid_cumsums <- lapply(cumsums, function(x) x[x <= n_clique])
  indices <- lapply(valid_cumsums, function(x) 1:length(x))
  set <- unlist(mapply("[", init_list, indices)) # set of numbers that can add up to n_clique
  set <- set[!is.na(set)] # some cleanup
  set
}

## get a random (match for a) combination of fitting cliques rather than the first, which is returned by match()
random_match <- function(combination, frequencies) {
  frequencies <- c(frequencies) # comes in as table, not vector
  df <- data.frame(frequencies)
  df$original_order <- 1:nrow(df)
  df <- df[sample(1:nrow(df)), ]
  df$original_order[match_only_once(combination, df$frequencies)]
}

# match() will return the same index repeatedly, but we can't use that because it is "used" once found
# so match each value individually sequentially and remove indices afterwards.
match_only_once <- function(combination, frequencies) {
  matches <- rep(NA, length(combination))
  for (i in seq_along(combination)) {
    matches[i] <- match(combination[i], frequencies)
    frequencies[matches[i]] <- NA
  }
  matches
}
