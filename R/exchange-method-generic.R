
#' Solve anticlustering using the modified exchange method
#'
#' @param data the data -- a N x N dissimilarity matrix or a N x M
#'     table of item features
#' @param K The number of cluster or an initial cluster assignment
#' @param obj_function the objective function. Takes as first argument
#'     a cluster assignment and as second argument the data set `data`
#'     (`data` is matrix, no `data.frame`).
#' @param categories A vector representing preclustering/categorical constraints
#'
#' @return The anticluster assignment
#'
#' @noRd
#'
#'

exchange_method <- function(data, K, obj_function, categories) {

  clusters <- initialize_clusters(NROW(data), K, categories)
  N <- nrow(data)
  best_total <- obj_function(data, clusters)
  for (i in 1:N) {
    # cluster of current item
    exchange_partners <- get_exchange_partners(clusters, i, categories)
    ## Do not use this item if there are zero exchange partners
    if (length(exchange_partners) == 0) {
      next
    }
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, length(exchange_partners))
    for (j in seq_along(exchange_partners)) {
      ## Swap item i with all legal exchange partners and check out objective
      comparison_objectives[j] <- update_objective_generic(
        data,
        clusters,
        i,
        exchange_partners[j],
        obj_function
      )
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # Which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      # Swap the elements
      clusters <- cluster_swap(clusters, i, swap)
      # Update best solution
      best_total <- best_this_round
    }
  }
  clusters
}

# Update the objective value - generic version taking any objective function
# param data: the data (distances/features)
# param clusters: a vector of clusters
# param i, j: the to be swapped items (indexes)
# param obj_function: The function that computes the objective value
# return: The objective after a hypothetical swap
update_objective_generic <- function(data, clusters, i, j, obj_function) {
  tmp_clusters <- cluster_swap(clusters, i, j)
  obj_function(data, tmp_clusters)
}

# Swap two items and return new arranged clusters
#
# param clusters: a vector of clusters
# param i, j: the to be swapped items (indexes of the items)
# return: the cluster vector after swapping items i and j
cluster_swap <- function(clusters, i, j) {
  group_i <- clusters[i]
  clusters[i] <- clusters[j]
  clusters[j] <- group_i
  clusters
}

# For a given item i, get all feasible exchange partners
#
# param clusters: a vector of clusters
# param i: the to be swapped item (index of the item)
# param categories: A vector of categories inducing constraints wrt
#                   which items are feasible exchange partners
# return: The indexes of the feasible exchange partners
get_exchange_partners <- function(clusters, i, categories) {
  N <- length(clusters)
  # are there categorical variables?
  if (!is.null(categories)) {
    # only exchange within the same group
    allowed_category <- categories == categories[i]
  } else {
    allowed_category <- rep(TRUE, N) # no constraint
  }

  ## Feasible exchange partners:
  # (a) are in different anticluster
  # (b) are in the same category/precluster
  exchange_partners <- (clusters != clusters[i]) & allowed_category
  (1:N)[exchange_partners]
}


# Initialize a cluster assignment.
# The initial assignment is either:
# (a) already passed via the K argument
# (b) a completely random assignment
# (c) a random assignment that balances the categories across anticlusters
#     (if categories are preclusters: each preclustered item is assigned
#      to a different anticlusters)
# param N: The number of cases
# param K: The number of clusters or an initial clustering assignment
# param categories: A vector of categories
# return: The initialized clusters
initialize_clusters <- function(N, K, categories) {
  if (length(K) > 1) {
    K <- to_numeric(K)
    return(K) # K is already an anticluster assignment
  }
  ## Initial assignment based on categorical constraints
  ## (categorical constraints may be preclustering constraints)
  if (argument_exists(categories)) {
    return(categorical_sampling(categories, K))
  }
  ## Initial random assignment unrestricted:
  sample(rep_len(1:K, length.out = N))
}
