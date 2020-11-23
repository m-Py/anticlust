
#' Wrapper for balanced clustering function in C
#' 
#' @param data N x M feature matrix
#' @param K The number of clusters
#' 
#' @return Clustering vector
#'
#' @noRd

foo_clust <- function(data, K) {
  data <- as.matrix(data)
  N <- nrow(data)
  if (N %% K != 0) {
    stop("NO! K BAD")
  }
  M <- ncol(data)
  order <- order(distances_from_centroid(data)) - 1 # used as index in C, so starts at 0
  results <- .C(
    "c_balanced_clustering", 
    as.double(data),
    as.integer(N),
    as.integer(M),
    as.integer(K),
    as.integer(order),
    vector = as.integer(rep(0, N)),
    mem_error = as.integer(0),
    PACKAGE = "anticlust"
  )
  if (results[["mem_error"]] == 1) {
    stop("Could not allocate enough memory.")
  }
  order_cluster_vector(results[["vector"]] + 1)
}
