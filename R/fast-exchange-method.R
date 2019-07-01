
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

get_neigbours <- function(features, k_neighbours, categories) {
  ncategories <- length(unique(categories))
  nns <- by(features, categories, RANN::nn2, k = k_neighbours + 1)
  nns <- lapply(nns, function(x) x$nn.idx)
  ## get original indices for each category
  ## (converts indices per category in global indices)
  get_index <- function(i, x) {
    which(categories == i)[x[[i]]]
  }
  idx <- lapply(1:ncategories, get_index, x = nns)
  dims <- lapply(nns, function(x) dim(x))
  ## change dimensions to original nearest neighbour structure
  adjust_dims <- function(a, b) {
    dim(a) <- dim(b)
    a
  }
  idx <- mapply(adjust_dims, idx, nns)
  idx <- do.call(rbind, idx)
  sort_by_col(idx, 1)
}
