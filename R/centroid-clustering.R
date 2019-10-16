#' Heuristic preclustering
#'
#' Computes equal-sized clusters based on the function
#' `centroid_anticlustering` by m.eik
#'
#' @param data N x M matrix of features or N x N distance matrix.
#' @param param N The sample size
#' @param param K the number of clusters
#' @return the cluster assignments as a vector
#'
#' @noRd
#'
#' @details The data argument can be distinguished by its class; either
#'     "features" or "distances". (both are an R matrix)
#'

centroid_clustering <- function(data, K) {
  N <- nrow(data)
  K <- N / K
  anticlusters <- centroid_anticlustering(data, k = K, as_vector = FALSE)
  mode(anticlusters) <- "numeric"
  ## sometimes the last item is not assigned, fix this:
  anticlusters[is.na(anticlusters)] <- which(!(1:N %in% anticlusters))
  preclusters <- cbind(c(anticlusters), rep(1:nrow(anticlusters), K))
  # return preclusters in correct order
  order_cluster_vector(preclusters[, 2][order(preclusters[, 1])])
}
