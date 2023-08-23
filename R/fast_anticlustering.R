
#' Fast anticlustering
#'
#' Increasing the speed of (k-means / k-plus) anticlustering by
#' conducting fewer exchanges during the optimization. This function
#' uses an adjusted exchange method in which the number of exchange
#' partners can be specified. Using fewer exchange partners can make
#' anticlustering applicable to quite large data sets.
#'
#' @param x A numeric vector, matrix or data.frame of data points.
#'     Rows correspond to elements and columns correspond to
#'     features. A vector represents a single numeric feature.
#' @param K How many anticlusters should be created. Alternatively:
#'     (a) A vector describing the size of each group, or (b) a vector
#'     of length \code{nrow(x)} describing how elements are assigned
#'     to anticlusters before the optimization starts.
#' @param k_neighbours The number of neighbours that serve as exchange
#'     partner for each element. Defaults to \code{Inf}, implying that
#'     each element is exchanged with each element in other groups.
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical constraints.
#' @param Optional argument. A list of length
#'     \code{NROW(x)} specifying for each element the indices of the
#'     elements that serve as exchange partners. If used, this
#'     argument overrides the \code{k_neighbours} argument. See
#'     examples.
#' @param backend Either "C" or "R", to use a C or R implementation of
#'     the anticlustering optimization algorithm. Since
#'     \code{anticlust} version 0.7.1, the faster C implementation is
#'     the default.
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
#' Papenberg, M. (2023). K-plus Anticlustering: An Improved k-means Criterion 
#' for Maximizing Between-Group Similarity. British Journal of Mathematical 
#' and Statistical Psychology. Advance online publication. 
#' https://doi.org/10.1111/bmsp.12315
#'
#'
#' @details
#'
#' This function was created to make anticlustering applicable to
#' large data sets (e.g., several 100,000 elements). It optimizes the
#' k-means objective because computing all pairwise distances as is
#' done when optimizing the "diversity" (i.e., the default in
#' \code{\link{anticlustering}}) is not feasible for very large data
#' sets (for about N > 20000 on my personal computer). Extending
#' \code{fast_anticlustering} to k-plus anticlustering is rather straight
#' forward by applying \code{\link{kplus_moment_variables}} on the input (see
#' examples).
#' 
#' The function \code{fast_anticlustering} employs a speed-optimized
#' exchange method, which is basically equivalent to \code{method =
#' "exchange"} in \code{\link{anticlustering}}, but (possibly) reduces
#' the number of exchanges that are investigated for each input
#' element. For each element, the potential exchange partners are
#' generated using a nearest neighbour search with the function
#' \code{\link[RANN]{nn2}} from the \code{RANN} package. Only the
#' nearest neighbours then serve as exchange partners. The number of
#' exchange partners per element has to be set using the argument
#' \code{k_neighbours}; by default, it is set to \code{Inf}, meaning
#' that all possible swaps are tested. This default must be changed by
#' the user for large data sets. If the default is not changed, you
#' can also just use the function \code{\link{anticlustering}}. More
#' exchange partners generally improve the quality of the results, but
#' also increase run time. Note that for very large data sets,
#' anticlustering generally becomes "easier" (even a random split may
#' yield satisfactory results), so using few exchange partners is
#' usually not a problem.
#' 
#' For a fixed number of exchange partners (specified using the
#' argument \code{k_neighbours}) the approximate run time of
#' \code{fast_anticlustering} is in O(M N^2), where N is the total
#' number of elements and M is the number of variables. The algorithm
#' \code{method = "exchange"} in \code{\link{anticlustering}} has a
#' run time of O(M N^3) because for each element, all other elements
#' serve as exchange partners. Thus, \code{fast_anticlustering} can
#' improve the run time by an order of magnitude as compared to the
#' standard exchange algorithm. The nearest neighbour search, which is
#' done in the beginning, only has O(N log(N)) run time and therefore
#' does not strongly contribute to the overall run time. It is
#' possible to suppress the nearest neighbour search by passing custom
#' exchange partners using the \code{exchange_partners} argument. The
#' examples show how to generate random exchange partners. When using
#' \code{exchange_partners}, it is not necessary that each element has
#' the same number of exchange partners (this is why the argument
#' \code{exchange_partners} has to be a list instead of matrix / data
#' frame).
#'
#' When setting the \code{categories} argument, exchange partners
#' (i.e., nearest neighbours) will be generated from the same
#' category. Note that when \code{categories} has multiple columns
#' (i.e., each element is assigned to multiple columns), each
#' combination of categories is treated as a distinct category by the
#' exchange method. You can also use
#' \code{\link{categories_to_binary}} to improve results for several
#' categorical variables, instead of using the argument
#' \code{categories}.
#'
#' @examples
#'
#' features <- iris[, - 5]
#'
#' ac_exchange <- fast_anticlustering(features, K = 3)
#'
#' ## The following call is (quasi) equivalent to the call above:
#' ac_exchange <- anticlustering(features, K = 3, objective = "variance")
#'
#' ## Improve run time by using fewer exchange partners:
#' ac_fast <- fast_anticlustering(features, K = 3, k_neighbours = 10)
#'
#' # Applying k-plus anticlustering with this function is straight forward,
#' # just use kplus_moment_variables() on the numeric input:
#' kplus_anticlusters <- fast_anticlustering(kplus_moment_variables(features, T = 2), K = 3)
#' mean_sd_tab(features, kplus_anticlusters) # Means and SDs are similar
#' 
#' # Working on several 1000 elements is very fast (Here n = 5000)
#' data <- matrix(rnorm(5000 * 2), ncol = 2)
#' groups <- fast_anticlustering(data, K = 2, k_neighbours = 2)
#' mean_sd_tab(data, groups)
#' 
#' # Use custom exchange partners, here: 10 random exchange partners for each element
#' n_exchange_partners <- 10
#' K <- 10
#' init <- sample(rep_len(1:K, nrow(features)))
#' groups_rnd_partners <- fast_anticlustering(
#'   features, 
#'   K = init, 
#'   exchange_partners = generate_exchange_partners(
#'     n_exchange_partners, 
#'     features = features, method = "random"
#'   )
#' )
#' 
#' # compare with using nearest neighbours as exchange partners (i.e., the default)
#' groups_nn_partners <- fast_anticlustering(features, K = init, k_neighbours = n_exchange_partners)
#' groups_all_partners <- fast_anticlustering(features, K = init)
#' 
#' variance_objective(features, groups_nn_partners)
#' variance_objective(features, groups_rnd_partners)
#' variance_objective(features, groups_all_partners)
#'

fast_anticlustering <- function(x, K, k_neighbours = Inf, categories = NULL, 
                                exchange_partners = NULL, backend = "C") {
  input_validation_anticlustering(x, K, "variance",
                                "exchange", FALSE, categories, NULL)
  validate_input(backend, "backend", objmode = "character", len = 1,
                 not_na = TRUE, not_function = TRUE, input_set = c("R", "C"))
  categories <- merge_into_one_variable(categories)
  if (!isTRUE(k_neighbours == Inf)) {
    validate_input(k_neighbours, "k_neighbours", objmode = "numeric", len = 1,
                   must_be_integer = TRUE, greater_than = 0, not_na = TRUE)
  }
  x <- as.matrix(x)
  N <- nrow(x)
  if (argument_exists(exchange_partners)) {
    validate_exchange_partners(exchange_partners, categories, N)
  } else {
    exchange_partners <- all_exchange_partners(x, k_neighbours, categories)
  }
  init <- initialize_clusters(N, K, categories)
  if (backend == "C") {
    # convert list of exchange partners to matrix for C; 
    # when doing that, possibly "fill" empty entries with "N+1" (in C -> N), if some
    # list elements are shorter than others (if categorical variables have been used and are
    # unevenly distributed)
    max_exchanges_partners <- max(lengths(exchange_partners))
    exchange_partners <- lapply(exchange_partners, function(x) c(x[1:length(x)], rep(N+1, max(0, max_exchanges_partners - length(x)))))
    # remove potential NAs
    exchange_partners <- lapply(exchange_partners, function(x) x[!is.na(x)]) 
    exchange_partners <- unname(t(t(as.data.frame(exchange_partners))))
    # `exchange_partners` is passed as matrix to C; there it is converted to a 1-dimensional "vector".
    # Here we pass it as a matrix where elements = columns; cols = exchange partners.
    # Ususally it would be the other way around, e.g. RANN returns the standard "N x k_neighbours"
    # format, but I pass a "k_neighbours x N" matrix to C, which is more easily handled there.
    return(c_anticlustering(x, init, categories = NULL, objective = "fast-kmeans", exchange_partners - 1))
  }
  fast_exchange_(x, init, exchange_partners)
}

validate_exchange_partners <- function(exchange_partners, categories, N) {
  if (!inherits(exchange_partners, "list")) {
    stop("Argument `exchange_partners` must be a list.")
  }
  if (argument_exists(categories)) {
    stop("Cannot use arguments `exchange_partners` and `categories` at the same time.")
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

#' Get exchange partners for fast_anticlustering()
#' 
#' @param n_exchange_partners The number of exchange partners per element
#' @param features The features for which nearest neighbours are sought. May be NULL if random exchange partners are generated.
#' @param N The number of elements for which exchange partners; can be \code{NULL} if \code{features} is passed (it is ignored if \code{features} is passed).
#' @param method Currently supports "RANN" and "random" 
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical constraints.
#'
#' @return A list of length \code{N}. Is usually used as input to the argument \code{exchange_partners} in \code{\link{fast_anticlustering}}.
#'   Then, the i'th element of the list contains the indices of the exchange partners that are used for the i'th element.
#' 
#' @export
#' 
#' @examples
#' 
#' generate_exchange_partners(5, N = 10, method = "random")
#' # may return less than 5 exchange partners if there are not enough members 
#' # of the same category: 
#' generate_exchange_partners(
#'   5, N = 10, 
#'   method = "random", 
#'   categories = cbind(schaper2019$room, schaper2019$frequency)
#' )
#' # using nearest neighbour search (unlike RANN::nn2, this does not 
#' # return the ID of the element itself as neighbour)
#' generate_exchange_partners(5, features = schaper2019[, 3:5], method = "RANN")[1:3]
#' # compare with RANN directly:
#' RANN::nn2(schaper2019[, 3:5], k = 5)$nn.idx[1:3, ]
#' 
generate_exchange_partners <- function(n_exchange_partners, features = NULL, N = NULL, method = "RANN", categories = NULL) {
  if (argument_exists(features)) {
    N <- NROW(features)
  }
  categories <- merge_into_one_variable(categories)
  validate_input(n_exchange_partners, "n_exchange_partners", objmode = "numeric", len = 1,
                 must_be_integer = TRUE, greater_than = 0, not_na = TRUE)
  if (argument_exists(features) && method == "RANN") {
    return(nearest_neighbours(features, n_exchange_partners, categories))
  } else if (method == "random") {
    init_list <- NULL
    if (argument_exists(categories)) {
      init_list <- list_idx_by_category(categories)
    }
    return(random_exchange_partners(N, n_exchange_partners, init_list))
  } else {
    stop("Argument method must be 'RANN' or 'random'. When method is 'RANN', the argument 'features' must be used.")
  }
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


# generate list of random exchange partners; will not generate itself as exchange partner
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


#' # remove the index of the element itself, so we ensure that each element has 
#' # 10 exchange partners excluding itself
