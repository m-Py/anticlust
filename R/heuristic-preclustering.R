#' Heuristic preclustering
#'
#' Computes equal-sized clusters based on the function
#' `centroid_anticlustering` by m.eik
#'
#' @param features N x M matrix of features
#' @param distances N x N distance matrix
#' @param param K the number of groups
#' @return the cluster assignments as a vector
#'
#' @noRd

centroid_preclustering <- function(data, distances, K) {
  if (argument_exists(data)) {
    N <- nrow(data)
  } else if (argument_exists(distances)) {
    distances <- as.matrix(distances)
    N <- nrow(distances)
  }
  K <- N / K
  anticlusters <- centroid_anticlustering(data, distances, k = K, as_vector = FALSE)
  mode(anticlusters) <- "numeric"
  ## sometimes the last item is not assigned, fix this:
  anticlusters[is.na(anticlusters)] <- which(!(1:N %in% anticlusters))
  preclusters <- cbind(c(anticlusters), rep(1:nrow(anticlusters), K))
  # return preclusters in correct order
  preclusters[, 2][order(preclusters[, 1])]
}

