
#' Fast anticlustering
#'
#' The most efficient way to solve anticlustering optimizing the
#' k-means variance criterion. Can be used for very large data sets.
#'
#' @param features A numeric vector, matrix or data.frame of data
#'     points.  Rows correspond to elements and columns correspond to
#'     features. A vector represents a single numeric feature.
#' @param K How many anticlusters should be created.
#' @param k_neighbours The number of neighbours that serve as exchange
#'     partner for each element. Defaults to Inf, i.e., each element
#'     is exchanged with each element in other groups.
#' @param categories A vector, data.frame or matrix representing one or
#'     several categorical constraints.
#'
#' @importFrom RANN nn2
#'
#' @export
#'
#' @examples
#'
#' start <- Sys.time()
#' features <- iris[, - 5]
#' ac_exchange <- fast_anticlustering(features, K = 3)
#' Sys.time() - start
#' by(features, ac_exchange, function(x) round(colMeans(x, na.rm = TRUE), 2))
#'
#' start <- Sys.time()
#' features <- schaper2019[, 3:6]
#' ac_exchange <- fast_anticlustering(features, K = 3, categories = schaper2019$room)
#' Sys.time() - start
#' by(features, ac_exchange, function(x) round(colMeans(x, na.rm = TRUE), 2))
#'

fast_anticlustering <- function(features, K, k_neighbours = Inf, categories = NULL) {
  anticlustering_(
    features,
    K = K,
    k_neighbours = k_neighbours,
    method = "fast-exchange",
    categories = categories
  )
}

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

fast_exchange_ <- function(data, clusters, categories, nearest_neighbors) {
  N <- nrow(data)
  best_total <- variance_objective_(clusters, data)
  centers <- cluster_centers(data, clusters)
  ## frequencies of each cluster are required for updating cluster centers:
  tab <- c(table(clusters))
  for (i in 1:N) {
    # cluster of current item
    cluster_i <- clusters[i]
    # are there categorical variables?
    exchange_partners <- nearest_neighbors[i, -1]
    if (!is.null(categories)) {
      exchange_partners <- exchange_partners[categories[exchange_partners] == categories[i]]
    }
    ## Do not change with other elements that are in the same cluster
    exchange_partners <- exchange_partners[clusters[exchange_partners] != cluster_i]
    ## Sometimes an exchange cannot take place
    if (length(exchange_partners) == 0) {
      next
    }
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, length(exchange_partners))
    for (j in seq_along(exchange_partners)) {
      ## Swap item i with all legal exchange partners and check out objective
      # (a) Determine clusters of to-be-swapped elements
      tmp_clusters <- clusters
      tmp_swap <- exchange_partners[j]
      cluster_j <- tmp_clusters[tmp_swap]
      tmp_clusters[i] <- cluster_j
      tmp_clusters[tmp_swap] <- cluster_i
      # (b) Update centers
      tmp_centers <- update_centers(centers, data, i, tmp_swap, cluster_i, cluster_j, tab)
      # (c) Compute distance from updated centers
      distances <- dist_from_centers(data, tmp_centers, squared = TRUE)
      # (d) Use two-column matrix to select distances that enter the
      #     objective function
      distances <- distances[cbind(1:nrow(distances), tmp_clusters)]
      # (e) Compute objective after exchange
      comparison_objectives[j] <- sum(distances)
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      ## Update centers
      centers <- update_centers(centers, data, i, swap, clusters[i], clusters[swap], tab)
      # swap the elements - update clusters
      clusters[i] <- clusters[swap]
      clusters[swap] <- cluster_i
      # update best solution
      best_total <- best_this_round
    }
  }
  clusters
}



#' Update a cluster center after swapping two elements
#'
#' @param centers The current cluster centers
#' @param features The features
#' @param i the index of the first element to be swapped
#' @param j the index of the second element to be swapped
#' @param cluster_i the cluster of element i
#' @param cluster_j the cluster of element j
#' @param tab A table of the cluster frequencies
#'
#' @details
#'
#' This should make the fast exchange method much faster, because
#' most time is spent on finding the cluster centers. After swapping
#' only two elements, it should be possible to update the two centers
#' very fast
#' @noRd

update_centers <- function(centers, features, i, j, cluster_i, cluster_j, tab) {
  ## First cluster: item i is removed, item j is added
  centers[cluster_i, ] <- centers[cluster_i, ] - (features[i, ] / tab[cluster_i]) + (features[j, ] / tab[cluster_j])
  ## Other cluster: item j is removed, item i is added
  centers[cluster_j, ] <- centers[cluster_j, ] + (features[i, ] / tab[cluster_i]) - (features[j, ] / tab[cluster_j])
  centers
}

#' Get neigbours for fast preclustering (by category)
#' @noRd
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

get_neighbours <- function(features, k_neighbours, categories) {

  if (argument_exists(categories)) {
    min_n <- min(table(categories))
    k_neighbours <- min(k_neighbours + 1, min_n)
  } else {
    k_neighbours <- min(nrow(features), k_neighbours + 1)
  }

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
  if (!argument_exists(categories)) {
    idx <- RANN::nn2(features, k = k_neighbours)$nn.idx
  } else {
    categories <- as.numeric(as.factor(categories))
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
    for (i in 1:length(idx)) {
      dim(idx[[i]]) <- dim(nns[[i]])
    }
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
    not_complete_idx <- cbind(not_complete, matrix(NA_partners, ncol = k_neighbours - 1))
    idx <- rbind(idx, not_complete_idx)
    colnames(idx) <- NULL
  }
  sort_by_col(idx, 1)
}
