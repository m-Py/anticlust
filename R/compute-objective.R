
#' Compute objective value for the (anti)clustering problem
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used in (anti)clustering.
#' @param anticlusters A vector representing the anticluster affiliation
#' @param objective The objective to be maximized, either "distance" or
#'     "variance".
#'
#' @return Scalar, the objective value; either the total
#'   within-cluster variance (if objective == "variance") or
#'   the total sum of within-cluster pointwise distances
#'   (if objective == "distance")
#'
#' @export
#'
#' @examples
#'
#' features <- matrix(rnorm(1000, 100, 15), ncol = 2)
#' n_elements <- 15
#' m_features <- 2
#' n_clusters <- 3
#' features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
#' clusters_exact <- clustering(features, n_clusters, method = "exact", standardize = FALSE)
#' clusters_heuristic <- clustering(features, n_clusters, method = "heuristic", standardize = FALSE)
#' # Exact clustering minimizes distance criterion
#' get_objective(features, clusters_exact, "distance")
#' get_objective(features, clusters_heuristic, "distance")
#' get_objective(features, clusters_exact, "variance")
#' # Heuristic clustering tries to minimize variance criterion
#' get_objective(features, clusters_heuristic, "variance")
#'
#' # Plot the clustering
#' par(mfrow = c(1, 2))
#' plot_clusters(features, clusters_exact)
#' plot_clusters(features, clusters_heuristic)
#'
#' @references
#' H. Späth, “Anticlustering: Maximizing the variance criterion,”
#' Control and Cybernetics, vol. 15, no. 2, pp. 213–218, 1986.
#'

get_objective <- function(features, anticlusters, objective) {
  features <- as.matrix(features)
  if (!objective %in% c("distance", "variance"))
    stop("Argument objective must be 'distance' or 'variance'.")
  if (objective == "distance")
    return(obj_value_distance(features, anticlusters))
  if (objective == "variance")
    return(obj_value_variance(features, anticlusters))
}


# Compute objective value for variance criterion
#
# @param features A data.frame, matrix or vector representing the
#     features that are used in the assignment.
# @param anticlusters A vector representing the anticluster affiliation
#
# @return Scalar, the total within-cluster variance.
#
#
obj_value_variance <- function(features, anticlusters) {
  ## 1. Compute cluster centers
  centers <- cluster_centers(features, anticlusters)
  ## 2. For each item, compute distance to each cluster center
  distances <- dist_from_centers(features, centers, squared = TRUE)
  ## 3. Use two-column matrix to select relevant distances
  distances <- distances[cbind(1:nrow(distances), anticlusters)]
  return(sum(distances))
}

# Objective value for the distance criterion
#
# @param features A data.frame, matrix or vector representing the
#     features that are used in the assignment.
# @param anticlusters A vector representing the anticluster affiliation
#
# @return Scalar, the total sum of within-cluster pointwise distances.
#
#
obj_value_distance <- function(features, anticlusters) {
  ## determine distances within each group
  distances <- by(features, anticlusters, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}

## Compute distance objective based on pre-computed distances
## (is better for complete enumeration and random sampling than `obj_value_distance`)
distance_objective <- function(distances, anticlusters, K) {
  sums_within <- rep(NA, K)
  for (k in 1:K) {
    ## is there a better / faster / less wasteful way to create all
    ## connections than expand.grid? (I only want those where one column
    ## has smaller value than the other)
    elements <- which(anticlusters == k)
    selection <- as.matrix(expand.grid(elements, elements))
    selection <- selection[selection[,1] < selection[,2], ]
    sums_within[k] <- sum(distances[selection])
  }
  return(sum(sums_within))
}

