
#' Solve anticlustering using the fast exchange method
#'
#' @param data the data -- an N x M table of item features
#' @param clusters An initial cluster assignment
#' @param categories A vector representing categorical constraints
#' @param nearest_neighbors A matrix of nearest neighbors given by RANN::nn2
#'
#' @return The anticluster assignment
#'
#' @noRd
#'
#'

fast_exchange_ <- function(data, clusters, categories, nearest_neighbors) {
  N <- nrow(data)
  best_total <- variance_objective_(clusters, data)
  n_damns <- 0
  for (i in 1:N) {
    # cluster of current item
    group_i <- clusters[i]
    # are there categorical variables?
    exchange_partners <- nearest_neighbors[i, -1]
    if (!is.null(categories)) {
      exchange_partners <- exchange_partners[categories[exchange_partners] == categories[i]]
    }
    ## do not change with item in the same group
    exchange_partners <- exchange_partners[clusters[exchange_partners] != group_i]
    ## Sometimes an exchange cannot take place
    if (length(exchange_partners) == 0) {
      n_damns <- n_damns + 1
      next
    }
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, length(exchange_partners))
    for (j in seq_along(exchange_partners)) {
      ## Swap item i with all legal exchange partners and check out objective
      tmp_clusters <- clusters
      tmp_clusters[i] <- tmp_clusters[exchange_partners[j]]
      tmp_clusters[exchange_partners[j]] <- group_i
      comparison_objectives[j] <- variance_objective_(tmp_clusters, data)
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      # swap the elements
      clusters[i] <- clusters[swap]
      clusters[swap] <- group_i
      # update best solution
      best_total <- best_this_round
    }
  }
  print(n_damns)
  clusters
}


#' Get neigbours for fast preclustering (by category)
#' @roRd
#'
#' @details
#'
#' Computes the k nearest neighbours for each input element using
#' RANN::nn2. Deals with NA: elements with NA are treated as
#' neighbours among each other because RANN::nn2 does not deal with
#' NA. Computes neighbors within each category because exchange
#' partners are only sought among neighbors and members of the same
#' category (it would not make sense to compute nearest neighbors
#' across categories because this function is only used as preprocessing
#' for the exchange algorithm)
#'

get_neigbours <- function(features, k_neighbours, categories) {

  ## indicices of non NA values; RANN::nn2 does not deal with NA
  bool_complete <- complete.cases(features)
  complete <- which(bool_complete)
  not_complete <- which(!bool_complete)
  ## exclude NA elements from nearest neighbor computation
  features <- features[complete, ]
  if (argument_exists(categories)) {
    categories <- categories[complete]
  }

  ## Compute nearest neighbors; within categories if categories
  ## are available!
  k_neighbours <- min(nrow(features), k_neighbours + 1)
  if (!argument_exists(categories)) {
    idx <- RANN::nn2(features, k = k_neighbours)$nn.idx
  } else {
    ncategories <- length(unique(categories))
    nns <- by(features, categories, RANN::nn2, k = k_neighbours)
    nns <- lapply(nns, function(x) x$nn.idx)
    ## get original indices for each category
    ## (converts indices per category in global indices)
    ## i is a category; x is the NN data structure
    get_index <- function(i, x) {
      which(categories == i)[x[[i]]]
    }
    idx <- lapply(1:ncategories, get_index, x = nns)
    ## change dimensions to original nearest neighbour structure
    adjust_dims <- function(a, b) {
      dim(a) <- dim(b)
      a
    }
    idx <- mapply(adjust_dims, idx, nns)
    idx <- do.call(rbind, idx)
    idx <- sort_by_col(idx, 1)
  }

  ## now deal with NA:
  if (sum(not_complete) >= 1) {
    ## recover original indices (before NA exclusion)
    dims <- dim(idx)
    idx <- complete[idx]
    dim(idx) <- dims ## recover original structure
    ## make elements with NA exchange partners among each other:
    NA_partners <- sample(not_complete, size = (k_neighbours - 1) * length(not_complete), replace = TRUE)
    not_complete_idx <- cbind(not_complete, matrix(NA_partners, ncol = (k_neighbours - 1)))
    idx <- rbind(idx, not_complete_idx)
    colnames(idx) <- NULL
  }
  sort_by_col(idx, 1)
}
