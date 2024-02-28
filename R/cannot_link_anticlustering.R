

#' Use cannot-link constraints with anticlustering
#' 
#' @param x An N x M feature matrix or N x N distance matrix
#' @param init_clusters vector of initial clusters or matrix with N columns where each row is an initial partition
#'   (usually a structure returned by optimal_dispersion() in the slot $groups)
#' @param cannot_link A 2 column matrix containing the indices of elements
#'     that must not be assigned to the same anticluster. (e.g., a matrix returned by
#'     optimal_dispersion() in the slot $edges)
#' @param objective "diversity", "average-diversity", "variance" or "kplus" (not dispersion!)
#' @param method "local-maximum" or "exchange" (TODO: "brusco"!)
#' 
#' @note
#' This function uses the average diversity objective if some groups are unequaled-sized
#' 
#' @noRd
cannot_link_anticlustering <- function(x, init_clusters, cannot_link, objective, method) {
  
  cannot_link <- as.matrix(cannot_link)
  
  if (objective == "kplus") {
    x <- kplus_moment_variables(x, 2)
    objective <- "variance"
  }
  
  if (objective == "variance") {
    x <- convert_to_distances(x)^2
  } else if (objective == "diversity") {
    x <- convert_to_distances(x)
  } 
  frequencies <- table(init_clusters)
  if (any(frequencies) != frequencies[1]) {
    objective <- "average-diversity"
  } else {
    objective <- "diversity"
  }


  # set cannot-link distances to large negative value so they cannot be linked
  x[rbind(cannot_link, t(apply(cannot_link, 1, rev)))] <- -(sum(x) + 1)
  
  
  ## special case of only one init partition...
  init_clusters <- as.matrix(init_clusters)
  if (ncol(init_clusters) == 1) {
    K <- c(init_clusters)
    init_clusters <- NULL
  } else {
    K <- init_clusters[1, ]
    init_clusters <- init_clusters - 1 # -1 for C
  }

  c_anticlustering(
    x, 
    K = K, 
    categories = NULL, 
    objective = objective, 
    local_maximum = ifelse(method == "local-maximum", TRUE, FALSE),
    exchange_partners = NULL,
    init_partitions = init_clusters
  )
  
}
