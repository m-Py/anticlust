
#' C implementation of anticlustering
#' 
#' @param data An N x M data matrix or N x N distance matrix.
#' @param K The number of clusters or an initial assignment of the elements 
#'     to clusters.
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical constraints. 
#' 
#' @noRd
#' 
c_anticlustering <- function(data, K, categories = NULL, objective, exchange_partners = NULL, local_maximum = FALSE, init_partitions = NULL) {
  
  clusters <- initialize_clusters(NROW(data), K, categories)

  # Ensure that `clusters` is a vector of integers between 0 and K-1
  clusters <- to_numeric(clusters) - 1
  # As convenience arguments for C, also pass 
  # a. the number of clusters (`K`)
  # b. how often each cluster occurs (`frequencies`)
  # c. the number of rows (`N`) and columns (`M`) in `data`
  K <- length(unique(clusters))
  frequencies <- table(clusters)
  N <- NROW(data)
  M <- NCOL(data)

  # Case: initial partitions may be passed as a matrix (rows = partitions; cols = elements)
  if (!argument_exists(init_partitions)) {
    R <- 1
    use_init_partitions <- 0
    init_partitions <- 0
  } else {
    R <- nrow(init_partitions)
    if (min(init_partitions) > 0) {
      stop("You forgot 0 indexing in C again.")
    }
    use_init_partitions <- 1
  }

  if (argument_exists(categories)) {
    USE_CATEGORIES <- TRUE
    categories <- merge_into_one_variable(categories) - 1
    N_CATS <- length(unique(categories))
    CAT_frequencies <- table(categories)
  } else {
    USE_CATEGORIES <- FALSE
    categories <- 0
    N_CATS <- 0
    CAT_frequencies <- 0
  }
  
  # Call C implementation of anticlustering
  if (objective == "variance") {
    results <- .C(
      "kmeans_anticlustering", 
      as.double(data),
      as.integer(N),
      as.integer(M),
      as.integer(K),
      as.integer(frequencies),
      clusters = as.integer(clusters),
      as.integer(USE_CATEGORIES),
      as.integer(N_CATS),
      as.integer(CAT_frequencies),
      as.integer(categories),
      mem_error = as.integer(0),
      PACKAGE = "anticlust"
    )
  } else if (objective %in% c("diversity", "distance", "average-diversity")) {
    if (objective != "average-diversity") {
      frequencies <- rep_len(1, K)
    }

    results <- .C(
      "distance_anticlustering", 
      as.double(convert_to_distances(data)),
      as.integer(N),
      as.integer(K),
      as.integer(frequencies),
      clusters = as.integer(clusters),
      as.integer(USE_CATEGORIES),
      as.integer(N_CATS),
      as.integer(CAT_frequencies),
      as.integer(categories),
      as.integer(local_maximum),
      as.integer(R),
      as.integer(use_init_partitions),
      as.integer(t(init_partitions)),
      mem_error = as.integer(0),
      PACKAGE = "anticlust"
    )
  } else if (objective == "dispersion") {
    results <- .C(
      "dispersion_anticlustering", 
      as.double(convert_to_distances(data)),
      as.integer(N),
      as.integer(K),
      clusters = as.integer(clusters),
      as.integer(USE_CATEGORIES),
      as.integer(N_CATS),
      as.integer(CAT_frequencies),
      as.integer(categories),
      mem_error = as.integer(0),
      PACKAGE = "anticlust"
    )
  } else if (objective == "fast-kmeans") {
    results <- .C(
      "fast_kmeans_anticlustering",
      as.double(data),
      as.integer(N),
      as.integer(M),
      as.integer(K),
      as.integer(frequencies),
      clusters = as.integer(clusters),
      as.integer(exchange_partners),
      as.integer(nrow(exchange_partners)),
      PACKAGE = "anticlust"
    )
    results[["mem_error"]] <- 0
  }

  if (results[["mem_error"]] == 1) {
    stop("Could not allocate enough memory.")
  }
  results[["clusters"]] + 1
}
