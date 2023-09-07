
#' Fast anticlustering
#'
#' Increasing the speed of (k-means / k-plus) anticlustering by (1) 
#' conducting fewer exchanges during the optimization and (2) using an alternative
#' formulation of the k-means objective. Makes anticlustering applicable to 
#' quite large data sets.
#'
#' @param x A numeric vector, matrix or data.frame of data points.
#'     Rows correspond to elements and columns correspond to
#'     features. A vector represents a single numeric feature.
#' @param K How many anticlusters should be created. Alternatively:
#'     (a) A vector describing the size of each group, or (b) a vector
#'     of length \code{nrow(x)} describing how elements are assigned
#'     to anticlusters before the optimization starts.
#' @param k_neighbours The number of nearest neighbours that serve as
#'     exchange partner for each element. See details.
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical constraints.
#' @param exchange_partners Optional argument. A list of length
#'     \code{NROW(x)} specifying for each element the indices of the
#'     elements that serve as exchange partners. If used, this
#'     argument overrides the \code{k_neighbours} argument. See
#'     examples.
#'
#' @importFrom RANN nn2
#'
#' @seealso
#'
#' \code{\link{anticlustering}}
#' 
#' \code{\link{kplus_moment_variables}}
#' 
#' \code{\link{categories_to_binary}}
#'
#' \code{\link{variance_objective}}
#' 
#' \code{\link{generate_exchange_partners}}
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' 
#' @references
#' 
#' 
#' Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition 
#' data sets into equivalent parts. Psychological Methods, 26(2), 
#' 161–174. https://doi.org/10.1037/met0000301.
#' 
#' Papenberg, M. (2023). K-plus Anticlustering: An Improved k-means
#' Criterion for Maximizing Between-Group Similarity. British Journal
#' of Mathematical and Statistical Psychology. Advance online
#' publication.  https://doi.org/10.1111/bmsp.12315
#' 
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213-218.
#'
#'
#' @details
#'
#' This function was created to make anticlustering applicable to
#' large data sets (e.g., several 100,000 elements). It optimizes the
#' k-means objective because computing all pairwise distances as is
#' done when optimizing the "diversity" (i.e., the default in
#' \code{\link{anticlustering}}) is not feasible for very large data
#' sets (for about N > 20000 on my personal computer). Using
#' \code{fast_anticlustering} for k-plus anticlustering is also
#' possible by applying \code{\link{kplus_moment_variables}} on the
#' input (and possibly by using the argument \code{exchange_partners},
#' see Examples).
#' 
#' The function \code{fast_anticlustering} employs a speed-optimized
#' exchange method, which is basically equivalent to \code{method =
#' "exchange"} in \code{\link{anticlustering}}, but reduces the number
#' of exchanges that are investigated for each input element. For each
#' element, the potential exchange partners are generated using a
#' nearest neighbour search with the function \code{\link[RANN]{nn2}}
#' from the \code{RANN} package. Only the nearest neighbours then
#' serve as exchange partners. The number of exchange partners
#' per element has to be set using the argument \code{k_neighbours}; by
#' default, it is set to \code{Inf}, meaning that all possible swaps are
#' tested. This default must be changed by the user for large data sets.
#' More exchange partners can improve the quality of
#' the results, but also increase run time. Note that for very large
#' data sets, anticlustering generally becomes "easier" (even a random
#' split may yield satisfactory results), so using few exchange
#' partners is usually not a problem. 
#' 
#' It is possible to specify custom exchange partners using the
#' argument \code{exchange_partners} instead of relying on the default
#' nearest neighbour search.  When using \code{exchange_partners}, it
#' is not necessary that each element has the same number of exchange
#' partners; this is why the argument \code{exchange_partners} has to
#' be a \code{list} instead of a \code{matrix} or
#' \code{data.frame}. Exchange partners can for example be generated
#' by \code{\link{generate_exchange_partners}} (see Examples), but a
#' custom list may also be used. Note that categorical constraints
#' induced via \code{categories} may not be respected during the
#' optimization if the \code{exchange_partners} argument allows
#' exchanges between members of different categories, so care must be
#' taken when combining the arguments \code{exchange_partners} and
#' \code{categories}.
#' 
#' In \code{anticlustering(..., objective = "variance")}, the run time
#' of computing the k-means objective is in O(M N), where N is the
#' total number of elements and M is the number of variables. This is
#' because the variance is computed as the sum of squared distances
#' between all data points and their cluster centers.  The function
#' \code{fast_anticlustering} uses a different - but equivalent -
#' formulation of the k-means objective, where the re-computation of
#' the objective only depends on K and M, but no longer on N. In
#' particular, it minimizes the weighted sum of squared distances between
#' cluster centroids and the overall data centroid; the distances
#' between all individual data points and their cluster center are not
#' computed (Späth, 1986). Using the different objective formulation 
#' reduces the run time by an
#' order of magnitude and makes k-means anticlustering applicable to
#' very large data sets (even in the millions). For a fixed number of
#' exchange partners (specified using the argument
#' \code{k_neighbours}), the approximate run time of
#' \code{fast_anticlustering} is in O(M N K). The algorithm
#' \code{method = "exchange"} in \code{\link{anticlustering}} with
#' \code{objective = "variance"} has a run time of O(M N^3). 
#' Thus, \code{fast_anticlustering} can improve the run time
#' by two orders of magnitude as compared to the standard exchange
#' algorithm. The nearest neighbour search, which is done in the
#' beginning, only has O(N log(N)) run time and does not strongly
#' contribute to the overall run time (and it is extremely fast in
#' practice). It is nevertheless possible to suppress the nearest
#' neighbour search by using the \code{exchange_partners} argument.
#'
#' When setting the \code{categories} argument, exchange partners
#' (i.e., nearest neighbours) will be generated from the same
#' category. Note that when \code{categories} has multiple columns,
#' each combination of these categories is treated as a distinct
#' category by the exchange method. You can also use
#' \code{\link{categories_to_binary}} to potentially improve results
#' for several categorical variables, instead of using the argument
#' \code{categories}.
#'
#' @examples
#'
#' ## Use fewer or more exchange partners to adjust speed (vs. quality tradeoff)
#' features <- iris[, - 5]
#' N <- nrow(features)
#' init <- sample(rep_len(1:3, N)) # same starting point for all calls:
#' groups1 <- fast_anticlustering(features, K = init) # default: all exchanges
#' groups2 <- fast_anticlustering(features, K = init, k_neighbours = 20) 
#' groups3 <- fast_anticlustering(features, K = init, k_neighbours = 2)
#' 
#' variance_objective(features, groups1)
#' variance_objective(features, groups2)
#' variance_objective(features, groups3)
#'
#' # K-plus anticlustering is straight forward when sticking with the default
#' # for k_neighbours
#' kplus_anticlusters <- fast_anticlustering(
#'   kplus_moment_variables(features, T = 2), 
#'   K = 3
#' )
#' mean_sd_tab(features, kplus_anticlusters)
#' 
#' # Some care is needed when applying k-plus using with this function 
#' # while using a reduced number of exchange partners generated in the 
#' # nearest neighbour search. Then we:
#' # 1) Use kplus_moment_variables() on the numeric input
#' # 2) Generate custom exchange_partners because otherwise nearest 
#' #    neighbours are internally selected based on the extended k-plus 
#' #    variables returned by kplus_moment_variables() 
#' #    (which does not really make sense)
#' kplus_anticlusters <- fast_anticlustering(
#'   kplus_moment_variables(features, T = 2), 
#'   K = 3,
#'   exchange_partners = generate_exchange_partners(120, features = features, method = "RANN")
#'  )
#' mean_sd_tab(features, kplus_anticlusters)
#' # Or we use random exchange partners: 
#' kplus_anticlusters <- fast_anticlustering(
#'   kplus_moment_variables(features, T = 2), 
#'   K = 3,
#'   exchange_partners = generate_exchange_partners(120, N = nrow(features), method = "random")
#' )
#' mean_sd_tab(features, kplus_anticlusters)
#' 
#' 
#' # Working on several 1000 elements is very fast (Here n = 10000, m = 2)
#' data <- matrix(rnorm(10000 * 2), ncol = 2)
#' start <- Sys.time()
#' groups <- fast_anticlustering(data, K = 5, k_neighbours = 5)
#' Sys.time() - start 
#'

fast_anticlustering <- function(x, K, k_neighbours = Inf, categories = NULL, 
                                exchange_partners = NULL) {
  input_validation_anticlustering(
    x, K, "variance", "exchange", FALSE, categories, NULL
  )
  categories <- merge_into_one_variable(categories)
  x <- as.matrix(x)
  N <- nrow(x)
  
  if (argument_exists(exchange_partners)) {
    validate_exchange_partners(exchange_partners, N)
  } else {
    if (!isTRUE(k_neighbours == Inf)) {
      validate_input(k_neighbours, "k_neighbours", objmode = "numeric", len = 1,
                     must_be_integer = TRUE, greater_than = 0, not_na = TRUE)
    }
    exchange_partners <- all_exchange_partners(x, k_neighbours, categories)
  }
  exchange_partners <- cleanup_exchange_partners(exchange_partners, N)

  c_anticlustering(
    x, initialize_clusters(N, K, categories), 
    categories = NULL, objective = "fast-kmeans", 
    exchange_partners - 1
  )
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
    return(list_idx_by_category(categories))
  }
  # Case 2: Everyone is potential exchange partner
  rep(list(1:N), N) 
}

# convert list of exchange partners to matrix for C; 
# when doing that, possibly "fill" empty entries with "N+1" (in C -> N), if some
# list elements are shorter than others (if categorical variables have been used and are
# unevenly distributed)
cleanup_exchange_partners <- function(exchange_partners, N) {
  max_exchanges_partners <- max(lengths(exchange_partners))
  exchange_partners <- lapply(exchange_partners, function(x) c(x[1:length(x)], rep(N+1, max(0, max_exchanges_partners - length(x)))))
  # remove potential NAs
  exchange_partners <- lapply(exchange_partners, function(x) x[!is.na(x)]) 
  exchange_partners <- unname(t(t(as.data.frame(exchange_partners))))
  # `exchange_partners` is passed as matrix to C; there it is converted to a 1-dimensional "vector".
  # Here we pass it as a matrix where elements = columns; cols = exchange partners.
  # Usually it would be the other way around, e.g. RANN returns the standard "N x k_neighbours"
  # format, but I pass a "k_neighbours x N" matrix to C, which is more easily handled there.
  exchange_partners
}

validate_exchange_partners <- function(exchange_partners, N) {
  if (!inherits(exchange_partners, "list")) {
    stop("Argument `exchange_partners` must be a list.")
  }
  if (length(exchange_partners) != N) {
    stop("If used, argument `exchange_partners` must have the same length as you have data points.")
  }
  if (any(lengths(exchange_partners)) > N) {
    stop("Argument `exchange_partners` has problems: Some elements have more exchange partners than there are data points.")
  }
  unlisted <- unlist(exchange_partners)
  if (any(unlisted) > N) {
    stop("Argument `exchange_partners` has problems: The maximum index is larger than the number of elements.")
  }
  if (any(unlisted) < 1) {
    stop("Argument `exchange_partners` has problems: Some indices are lower than 1.")
  }
  if (any(as.integer(unlisted) != unlisted)) {
    stop("Argument `exchange_partners` has problems: Only integer indices are allowed.")
  }
}


#' Get exchange partners for fast_anticlustering()
#' 
#' @param n_exchange_partners The number of exchange partners per
#'     element
#' @param N The number of elements for which exchange partners; can be
#'     \code{NULL} if \code{features} is passed (it is ignored if
#'     \code{features} is passed).
#' @param features The features for which nearest neighbours are
#'     sought if \code{method = "RANN"}.  May be NULL if random
#'     exchange partners are generated.
#' @param method Currently supports "random" (default), "RANN" and
#'     "restricted_random". See details.
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical constraints.
#'
#' @return A list of length \code{N}. Is usually used as input to the
#'     argument \code{exchange_partners} in
#'     \code{\link{fast_anticlustering}}.  Then, the i'th element of
#'     the list contains the indices of the exchange partners that are
#'     used for the i'th element.
#' 
#' @export
#' 
#' @details
#' 
#' The \code{method = "RANN"} generates exchange partners using a
#' nearest neighbour search via \code{\link[RANN]{nn2}} from the
#' \code{RANN} package; \code{methode = "restricted_random"} generates
#' random exchange partners but ensures that for each element, no
#' duplicates are generated and that the element itself does not occur
#' as exchange partner (this is the slowest method, and I would not
#' recommend it for large N); \code{method = "random"} (default) does
#' not impose these restrictions and generates unrescricted random
#' partners (it may therefore generate duplicates and the element
#' itself as exchange partner).
#' 
#' When setting the \code{categories} argument and using \code{method
#' = "RANN"}, exchange partners (i.e., nearest neighbours) will be
#' generated from the same category; \code{methode =
#' "restricted_random"} will also adhere to categorical constraints
#' induced via \code{categories} (i.e. each element only receives
#' exchange partners from the same category as itself); \code{methode
#' = "random"} cannot incoorporate categorical restrictions.
#' 
#' 
#' @examples
#' 
#' # Restricted random method generates no duplicates per element and cannot return 
#' # the element itself as exchange partner
#' generate_exchange_partners(5, N = 10, method = "restricted_random")
#' # "random" simply randomizes with replacement and without restrictions
#' # (categorical restrictions are also not possible; is much faster for large data sets)
#' generate_exchange_partners(5, N = 10, method = "random")
#' # May return less than 5 exchange partners if there are not enough members 
#' # of the same category: 
#' generate_exchange_partners(
#'   5, N = 10, 
#'   method = "restricted_random", 
#'   categories = cbind(schaper2019$room, schaper2019$frequency)
#' )
#' # using nearest neighbour search (unlike RANN::nn2, this does not 
#' # return the ID of the element itself as neighbour)
#' generate_exchange_partners(5, features = schaper2019[, 3:5], method = "RANN")[1:3]
#' # compare with RANN directly:
#' RANN::nn2(schaper2019[, 3:5], k = 6)$nn.idx[1:3, ] # note k = 6
#' 
generate_exchange_partners <- function(n_exchange_partners, N = NULL, features = NULL, method = "random", categories = NULL) {
  if (argument_exists(features)) {
      validate_data_matrix(features)
      N <- NROW(features)
  }
  categories <- merge_into_one_variable(categories)
  validate_input(n_exchange_partners, "n_exchange_partners", objmode = "numeric", len = 1,
                 must_be_integer = TRUE, greater_than = 0, not_na = TRUE)
  validate_input(method, "method", objmode = "character", len = 1,
                 not_na = TRUE, not_function = TRUE, input_set = c("RANN", "random", "restricted_random"))
  if (argument_exists(features) && method == "RANN") {
    return(nearest_neighbours(features, n_exchange_partners, categories))
  } else if (method == "restricted_random") {
    init_list <- NULL
    if (argument_exists(categories)) {
      init_list <- list_idx_by_category(categories)
    }
    return(random_exchange_partners(N, n_exchange_partners, init_list))
  } else if (method == "random" && !argument_exists(categories)) {
    return(completely_random_exchange_partners(N, n_exchange_partners))
  }
  else {
    stop("Argument method must be 'RANN' or 'random' or 'restricted_random'.\n'",
         "When method is 'RANN', the argument 'features' must be used.\n'",
         "When method is 'random', you cannot use the argument 'categories'.")
  }
}

list_idx_by_category <- function(categories) {
  category_ids <- lapply(1:max(categories), function(i) which(categories == i))
  category_ids[categories]
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
  remove_self(idx)
}

# Convert a matrix to list - each row becomes list element
matrix_to_list <- function(x) {
  as.list(as.data.frame(t(x)))
}

# Merge a list of lists into one list
merge_lists <- function(list_of_lists) {
  do.call(c, list_of_lists)
}

completely_random_exchange_partners <- function(N, n_exchange_partners) {
  ret <- matrix_to_list(matrix(sample(N, size = N * n_exchange_partners, replace = TRUE), ncol = n_exchange_partners))
  unname(ret)
}

# generate list of random exchange partners; will not generate itself as exchange partner
# and has no duplicates for each element
random_exchange_partners <- function(N, n_exchange_partners, init_partners = NULL) {
  if (is.null(init_partners)) {
    init_partners <- rep(list(1:N), N) 
  }
  init_partners <- remove_self(init_partners)
  partners <- lapply(init_partners, function(x) sample_(x)[1:n_exchange_partners])
  # possibly remove NAs
  lapply(partners, function(x) x[!is.na(x)]) 
}

# remove in a list for each element the ID of the list elements position 
remove_self <- function(idx_list) {
  N <- length(idx_list)
  lapply(1:N, function(i) idx_list[[i]][idx_list[[i]] != i])
}




#### THIS IS AN OLD AN OBSOLETE R IMPLEMENTATION; I AM KEEPING IT FOR THE TEST FILES
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

############ END OLD OBSOLETE IMPLEMENTATION
