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
  anticlusters <- as_numeric_anticlusters(data, anticlusters)
  # for each element, determine the cluster (anticlusters is a vector
  # with a meaningful order, caused by the matrix that is returned
  # by `centroid_anticlustering`
  preclusters  <- cbind(anticlusters, rep(1:(N / K), K))
  # return preclusters in correct order
  order_cluster_vector(preclusters[, 2][order(preclusters[, 1])])
}


# Sometimes the last element is not assigned to a cluster, fix this
# through a postprocessing. Returns a numeric vector of indexes of
# elements. The order of the elements can corresponds to the clusters
# each element has been assigned to.
as_numeric_anticlusters <- function(data, anticlusters) {
  # centroid_anticlustering returns a matrix of row names (if the input
  # had rownames) or a character matrix of numeric indices. each row
  # in the matrix corresponds to a cluster.
  id_numeric <- 1:nrow(data)
  if (is.null(rownames(data))) {
    idx <- as.character(id_numeric)
  }
  else {
    idx <- rownames(data)
  }
  names(id_numeric) <- idx
  anticlusters[anticlusters == ""] <- idx[!(idx %in% anticlusters)]
  ## return numeric indices of anticlusters
  unname(id_numeric[c(anticlusters)])
}
