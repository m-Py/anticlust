

#' Use cannot-link constraints with anticlustering
#' 
#' @param x An N x M feature matrix or N x N distance matrix
#' @param init_clusters vector of initial clusters or matrix with N columns where each row is an initial partition
#'   (usually a structure returned by optimal_dispersion() in the slot $groups)
#' @param cannot_link A 2 column matrix containing the indices of elements
#'     that must not be assigned to the same anticluster. (e.g., a matrix returned by
#'     optimal_dispersion() in the slot $edges)
#' @param objective "diversity", "average-diversity", "variance" or "kplus" (not dispersion!)
#' @param method "local-maximum" or "exchange", "brusco"
#' 
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
  } else if (grepl("diversity", objective)) {
    x <- convert_to_distances(x)
  } 
  frequencies <- table(init_clusters)
  if (objective == "variance") {
    objective <- "average-diversity"
  }
  stopifnot(grepl("diversity", objective)) # must be diversity or average-diversity here

  ## special case of only one init partition...
  init_clusters <- as.matrix(init_clusters)
  if (ncol(init_clusters) == 1) {
    K <- c(init_clusters)
    init_clusters_bils <- K # unfortunately, these functions have different interfaces...
    init_clusters <- NULL
  } else {
    K <- init_clusters[1, ]
    init_clusters_bils <- init_clusters
    init_clusters <- init_clusters - 1 # -1 for C
  }
  
  if (method == "brusco") {
    return(BILS_CANNOT_LINK(x, init_clusters_bils, cannot_link, objective))
  }
  
  # For LCW: Set cannot-link distances to large negative value so they cannot be linked
  x[cleanup_cannot_link_indices(cannot_link)] <- -(sum(x) + 1)

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


# Wrapper for BILS method that preserves optimal dispersion and potentially uses multiple initial partitions
BILS_CANNOT_LINK <- function(distances, init_clusters, cannot_link, objective) {
  N <- nrow(distances)
  multiple_partitions_as_input <- is.matrix(init_clusters)
  cl_matrix <- matrix(1, ncol = N, nrow = N)
  cl_matrix[cleanup_cannot_link_indices(cannot_link)] <- -1
  bicriterion_anticlustering(
    distances, 
    K = if (multiple_partitions_as_input) init_clusters[1, ] else init_clusters,
    R = if (multiple_partitions_as_input) rep(ceiling(nrow(init_clusters)/2), 2) else c(1, 1),
    init_partitions = if (multiple_partitions_as_input) init_clusters[1:ceiling(nrow(init_clusters)/2),] else NULL,
    dispersion_distances = cl_matrix,
    average_diversity = objective == "average-diversity",
    return = "best-dispersion"
  )
}

cleanup_cannot_link_indices <- function(cannot_link) {
  cannot_link <- rbind(cannot_link, t(apply(cannot_link, 1, rev))) # use (1, 2) and (2, 1)
  cannot_link[!duplicated(cannot_link), ] # but do not use (1, 2), (2, 1), (1, 2) and (2, 1)
} 
