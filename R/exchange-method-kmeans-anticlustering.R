
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
                                "exchange", FALSE, categories, NULL)
  categories <- merge_into_one_variable(categories)
  if (!isTRUE(k_neighbours == Inf)) {
    validate_input(k_neighbours, "k_neighbours", objmode = "numeric", len = 1,
                   must_be_integer = TRUE, greater_than = 0, not_na = TRUE)
  }
  x <- as.matrix(x)
  exchange_partners <- all_exchange_partners(x, k_neighbours, categories)
  init <- initialize_clusters(nrow(x), K, categories)
  fast_exchange_(x, init, exchange_partners)
}

#' Solve anticlustering using the fast exchange method
#'
#' @param data the data -- an N x M table of item features
#' @param clusters An initial cluster assignment
#' @param all_exchange_partners A list of exchange partners
#'
#' @return The anticluster assignment
#'
#' @noRd
#'

fast_exchange_ <- function(data, clusters, all_exchange_partners) {
  N <- nrow(data)
  best_total <- variance_objective_(clusters, data)
  centers <- cluster_centers(data, clusters)
  distances <- dist_from_centers(data, centers, squared = TRUE)
  
  ## frequencies of each cluster are required for updating cluster centers:
  tab <- c(table(clusters))
  for (i in 1:N) {
    # cluster of current item
    cluster_i <- clusters[i]
    # get exchange partners for item i
    exchange_partners <- all_exchange_partners[[i]]
    # exchange partners are not in the same cluster:
    exchange_partners <- exchange_partners[clusters[exchange_partners] != clusters[i]]
    # Sometimes an exchange cannot take place
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
      # (b) Swap the elements
      tmp_clusters[i] <- cluster_j
      tmp_clusters[tmp_swap] <- cluster_i
      # (c) Update cluster centers after swap
      tmp_centers <- update_centers(centers, data, i, tmp_swap, cluster_i, cluster_j, tab)
      # (d) Update distances from centers after swap
      tmp_distances <- update_distances(data, tmp_centers, distances, cluster_i, cluster_j)
      # (e) Compute objective after swap
      comparison_objectives[j] <- sum(tmp_distances[cbind(1:nrow(tmp_distances), tmp_clusters)])
    }
    ## If an improvement of the objective occured, do the swap
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      # Update cluster centers
      centers <- update_centers(centers, data, i, swap, clusters[i], clusters[swap], tab)
      # Update distances
      distances <- update_distances(data, centers, distances, cluster_i, clusters[swap])
      # Actually swap the elements - i.e., update clusters
      clusters[i] <- clusters[swap]
      clusters[swap] <- cluster_i
      # Update best solution
      best_total <- best_this_round
    }
  }
  clusters
}

#' Recompute distances from cluster centers after swapping two elements
#' @param distances distances from cluster centers per element (old)
#' @param cluster_i the cluster of element i
#' @param cluster_j the cluster of element j
#' @return The new distances
#' @noRd
update_distances <- function(features, centers, distances, cluster_i, cluster_j) {
  for (k in c(cluster_i, cluster_j)) {
    distances[, k] <- colSums((t(features) - centers[k,])^2)
  }
  distances
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
  centers[cluster_i, ] <- centers[cluster_i, ] - (features[i, ] / tab[cluster_i]) + (features[j, ] / tab[cluster_i])
  ## Other cluster: item j is removed, item i is added
  centers[cluster_j, ] <- centers[cluster_j, ] + (features[i, ] / tab[cluster_j]) - (features[j, ] / tab[cluster_j])
  centers
}

#' Get exchange partners for k-means anticlustering 
#'
#' @details
#'
#' Computes the k nearest neighbours for each input element using
#' RANN::nn2. If no nearest neighbours are required, argument 
#' `k_neighbours` will be `Inf`. May compute nearest neighbors within
#' categories. 
#' 
#' @return A list of length `N`. Each element is a vector
#'    of exchange partners (that may be nearest neighbors).
#'
#' @noRd


all_exchange_partners <- function(features, k_neighbours, categories) {
  # Case 1: no nearest neighbor search needed
  if (is.infinite(k_neighbours)) { 
    return(all_exchange_partners_(nrow(features), categories))
  }
  # Case 2: NN search needed
  return(nearest_neighbours(features, k_neighbours, categories))
}

# Generate all possible exchange partners 
all_exchange_partners_ <- function(N, categories) {
  # Case 1: Exchange partners are from the same category
  if (argument_exists(categories)) {
    category_ids <- lapply(1:max(categories), function(i) which(categories == i))
    return(category_ids[categories])
  }
  # Case 2: Everyone is potential exchange partner
  rep(list(1:N), N) 
}

# Generate exchange partners via nearest neighbor search using RANN::nn2
nearest_neighbours <- function(features, k_neighbours, categories) {
  if (!argument_exists(categories)) {
    idx <- matrix_to_list(RANN::nn2(features, k = min(k_neighbours + 1, nrow(features)))$nn.idx)
  } else {
    # compute nearest neighbors within each category
    nns <- list()
    # track the indices when dividing by category. problem is:
    # in nearest neighbour search, there is no guarantee that 
    # the nearest neighbour is the element itself; if it were, 
    # restoring original order would be easy
    new_order <- list()
    for (i in 1:max(categories)) {
      tmp_indices <- which(categories == i)
      new_order[[i]] <- tmp_indices
      tmp_features <- features[tmp_indices, , drop = FALSE]
      tmp_nn <- RANN::nn2(tmp_features, k = min(k_neighbours + 1, nrow(tmp_features)))$nn.idx
      # get original index per category
      nns[[i]] <- which(categories == i)[tmp_nn]
      # restore matrix structure, gets lost 
      dim(nns[[i]]) <- dim(tmp_nn)
    }
    # per category, convert matrix to list
    idx_list <- lapply(nns, matrix_to_list)
    # in the end, merge all lists into 1 list
    idx <- merge_lists(idx_list)
    # restore original order after dividing by category
    original_order <- order(unlist(new_order))
    idx <- idx[original_order]
  }
  # return as list
  idx
}

# Convert a matrix to list - each row becomes list element
matrix_to_list <- function(x) {
  as.list(as.data.frame(t(x)))
}

# Merge a list of lists into one list
merge_lists <- function(list_of_lists) {
  do.call(c, list_of_lists)
}
