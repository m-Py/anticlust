#' Heuristic preclustering
#'
#' Computes equal-sized clusters based on the function
#' `centroid_anticlustering` by m.eik
#'
#' @param features N x M matrix of features
#' @param distances N x N distance matrix
#' @param param N The sample size
#' @param param K the number of clusters
#' @return the cluster assignments as a vector
#'
#' @noRd

centroid_clustering <- function(features, distances, N, K) {
  K <- N / K
  if (argument_exists(features)) {
    rownames(features) <- NULL
  }
  if (argument_exists(distances)) {
    rownames(distances) <- NULL
  }
  anticlusters <- centroid_anticlustering(features, distances, k = K, as_vector = FALSE)
  mode(anticlusters) <- "numeric"
  ## sometimes the last item is not assigned, fix this:
  anticlusters[is.na(anticlusters)] <- which(!(1:N %in% anticlusters))
  preclusters <- cbind(c(anticlusters), rep(1:nrow(anticlusters), K))
  # return preclusters in correct order
  preclusters[, 2][order(preclusters[, 1])]
}

