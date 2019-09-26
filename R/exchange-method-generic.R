
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

  clusters <- initialize_clusters(data, K, obj_function, categories)
  N <- nrow(data)
  best_total <- obj_function(clusters, data)
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
      comparison_objectives[j] <- update_objective_generic(data, clusters, i, exchange_partners[j], obj_function)
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
  obj_function(tmp_clusters, data)
}

# Swap two items and return new arranged clusters
#
# param clusters: a vector of clusters
# param i, j: the to be swapped items (indexes)
# return: the cluster vector after swapping items i and j
cluster_swap <- function(clusters, i, j) {
  group_i <- clusters[i]
  clusters[i] <- clusters[j]
  clusters[j] <- group_i
  clusters
}

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
  # (c) are NA (and in the same precluster/category)
  exchange_partners <- (clusters != clusters[i]) & allowed_category
  ## NA is converted to TRUE. This works because: An item only has NA
  ## if its cluster is currently NA and has the same category/precluster
  ## (case (d) from above)
  exchange_partners[is.na(exchange_partners)] <- TRUE
  (1:N)[exchange_partners]
}


initialize_clusters <- function(data, K, obj_function, categories) {
  if (length(K) > 1) {
    return(K) # K is already an anticluster assignment
  }
  ## Initial assignment based on categorical constraints
  ## (categorical constraints may be preclustering constraints)
  if (argument_exists(categories)) {
    return(random_sampling(data, K, obj_function, nrep = 1, categories))
  }
  ## Initial random assignment unrestricted:
  sample(rep_len(1:K, length.out = nrow(data)))
}