
#' Fast anticlustering
#'
#' The most efficient way to solve anticlustering optimizing the
#' k-means variance criterion with an exchange method. Can be used for
#' very large data sets.
#'
#' @param x A numeric vector, matrix or data.frame of data
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
#' @seealso
#'
#' \code{\link{anticlustering}}
#'
#' \code{\link{variance_objective}}
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @details
#'
#' This function was created to make anticlustering applicable
#' to large data sets (e.g., 100,000 elements). It optimizes the k-means
#' variance objective because computing all pairwise distances is not
#' feasible for many elements. Additionally, this function employs a
#' speed-optimized exchange method. For each element, the potential
#' exchange partners are generated using a nearest neighbor search with the
#' function \code{\link[RANN]{nn2}} from the \code{RANN} package. The nearest
#' neighbors then serve as exchange partners. This approach is inspired by the
#' preclustering heuristic according to which good solutions are found
#' when similar elements are in different sets---by swapping nearest
#' neighbors, this will often be the case. The number of exchange partners
#' per element has to be set using the argument \code{k_neighbours}; by
#' default, it is set to \code{Inf}, meaning that all possible swaps are
#' tested. This default must be changed by the user for large data sets.
#' More exchange partners generally improve the output, but also increase
#' run time.
#'
#' When setting the \code{categories} argument, exchange partners will
#' be generated from the same category. Note that when
#' \code{categories} has multiple columns (i.e., each element is
#' assigned to multiple columns), each combination of categories is
#' treated as a distinct category by the exchange method.
#'
#' @examples
#'
#'
#' features <- iris[, - 5]
#'
#' start <- Sys.time()
#' ac_exchange <- fast_anticlustering(features, K = 3)
#' Sys.time() - start
#'
#' ## The following call is equivalent to the call above:
#' start <- Sys.time()
#' ac_exchange <- anticlustering(features, K = 3, objective = "variance")
#' Sys.time() - start
#'
#' ## Improve run time by using fewer exchange partners:
#' start <- Sys.time()
#' ac_fast <- fast_anticlustering(features, K = 3, k_neighbours = 10)
#' Sys.time() - start
#'
#' by(features, ac_exchange, function(x) round(colMeans(x), 2))
#' by(features, ac_fast, function(x) round(colMeans(x), 2))
#'

fast_anticlustering <- function(x, K, k_neighbours = Inf, categories = NULL) {
  input_validation_anticlustering(x, K, "variance",
                                "exchange", FALSE, categories)

  if (!isTRUE(k_neighbours == Inf)) {
    validate_input(k_neighbours, "k_neighbours", objmode = "numeric", len = 1,
                   must_be_integer = TRUE, greater_than = 0, not_na = TRUE)
  }
  x <- as.matrix(x)
  neighbours <- get_neighbours(x, k_neighbours, categories)
  init <- initialize_clusters(nrow(x), K, categories)
  fast_exchange_(x, init, categories, neighbours)
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
    exchange_partners <- get_exchange_partners_kmeans(clusters, i, nearest_neighbors, categories)
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

get_exchange_partners_kmeans <- function(clusters, i, nearest_neighbors, categories) {
  exchange_partners <- nearest_neighbors[i, -1]
  if (!is.null(categories)) {
    exchange_partners <- exchange_partners[categories[exchange_partners] == categories[i]]
  }
  ## Do not change with other elements that are in the same cluster
  exchange_partners[clusters[exchange_partners] != clusters[i]]
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
#' @importFrom stats complete.cases
#'
#' @noRd


get_neighbours <- function(features, k_neighbours, categories) {

  if (argument_exists(categories)) {
    min_n <- min(table(categories))
    k_neighbours <- min(k_neighbours + 1, min_n)
  } else {
    k_neighbours <- min(nrow(features), k_neighbours + 1)
  }

  ## Compute nearest neighbors; within categories if categories
  ## are available!
  if (!argument_exists(categories)) {
    idx <- RANN::nn2(features, k = k_neighbours)$nn.idx
  } else {
    categories <- to_numeric(categories)
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

  sort_by_col(idx, 1)
}
