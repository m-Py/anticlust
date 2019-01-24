
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
#' n_anticlusters <- 4
#' # Precluster cases
#' n_preclusters <- nrow(features) / n_anticlusters
#' preclusters <- equal_sized_kmeans(features, n_preclusters)
#' # Use preclustering as resticting information in anticlustering
#' anticlusters <- heuristic_anticlustering(features, preclusters, objective = "distance")
#' get_objective(features, anticlusters, "distance")
#' get_objective(features, anticlusters, "variance")
#' anticlusters <- heuristic_anticlustering(features, preclusters, objective = "variance")
#' get_objective(features, anticlusters, "distance")
#' get_objective(features, anticlusters, "variance")
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


#' Compute objective value for variance criterion
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used in the assignment.
#' @param anticlusters A vector representing the anticluster affiliation
#'
#' @return Scalar, the total within-cluster variance.
#'
#'
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
#' @return Scalar, the total sum of within-cluster pointwise distances.
#'
#'
obj_value_distance <- function(features, anticlusters) {
  ## determine distances within each group
  distances <- by(features, anticlusters, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}
