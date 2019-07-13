
#' Internal heuristic method for anticlustering
#'
#' @param preclusters A preclustering vector
#' @param k_neighbours The number of neighbours for the fast
#'     exchange method
#'
#' @noRd
#'
heuristic_anticlustering <- function(data, K, obj_function,
                                     method, preclusters, nrep,
                                     categories, parallelize,
                                     seed, k_neighbours) {
  if (method == "sampling" || method == "heuristic") {
    # For random sampling, we cannot apply preclustering and categorical
    # constraints at the same time
    if (argument_exists(categories)) {
      preclusters <- NULL
    }
    return(random_sampling(data, K, preclusters,
                           obj_function, nrep, categories,
                           parallelize, seed,
                           ncores = NULL))
  } else if (method == "exchange") {
    return(exchange_method(data, K, obj_function, categories, preclusters))
  } else if (method == "fast-exchange") {
    ## fast exchange only
    neighbours <- get_neighbours(data, k_neighbours, categories)
    clusters <- random_sampling(data, K, NULL, obj_function,
                                nrep = 1, categories, FALSE,
                                NULL, NULL)
    fast_exchange_(data, clusters, categories, neighbours)
  }
}
