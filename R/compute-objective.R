
#' Compute objective value for variance criterion
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used in the assignment.
#' @param anticlusters A vector representing the anticluster affiliation
#'
#' @return Scalar, the total within-cluster variance.
#'
#' @noRd

obj_value_variance <- function(features, anticlusters) {
  ## 1. Compute cluster centers
  centers <- cluster_centers(features, anticlusters)
  ## 2. For each item, compute distance to each cluster center
  distances <- dist_from_centers(features, centers, squared = TRUE)
  ## 3. Use two-column matrix to select relevant distances
  distances <- distances[cbind(1:nrow(distances), anticlusters)]
  return(sum(distances))
}

#' Objective value for the distance criterion
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used in the assignment.
#' @param anticlusters A vector representing the anticluster affiliation
#'
#' @return Scalar, the total sum of within-cluster distances (based
#'     on the Euclidean distance).
#'
#' @importFrom stats dist
#'
#' @noRd

obj_value_distance <- function(features, anticlusters) {
  ## determine distances within each group
  distances <- by(features, anticlusters, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}


## Distance objective based on pre-computed distances
## (is better for complete enumeration than `obj_value_distance`)
distance_objective <- function(distances, anticlusters, K) {
  sums_within <- rep(NA, K)
  for (k in 1:K) {
    elements <- which(anticlusters == k)
    selection <- unique_combinations(elements)
    sums_within[k] <- sum(distances[selection])
  }
  return(sum(sums_within))
}

## This variant to create unique combinations is *much* faster
## (especially for larger N) than using utils::combn even though it
## seems that a lot of unnecessary combinations are generated using
## expand.grid
unique_combinations <- function(elements) {
  selection <- as.matrix(expand.grid(elements, elements))
  selection[selection[, 1] < selection[, 2], ]
}
