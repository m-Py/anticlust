
#' Anticlustering based on a heuristic
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used.
#' @param clustering A vector representing the preclustering of
#'     elements. See Details.
#' @param objective The objective to be maximized, either
#'     "distance" (default) or "variance".
#' @param nrep The number of repetitions tried when assigning
#'     elements to anticlusters.
#'
#' @return A vector representing the anticlustering.
#'
#' @details The heuristic approach to anticlustering forbids elements
#'     belonging to the same precluster to be assigned to the same
#'     anticluster. The preclustering should be accomplished by one of the
#'     clustering functions, `equal_sized_cluster_editing` (an exact
#'     method that minimizes distance criterion under the restriction of
#'     equal group sizes) or `equal_sized_kmeans` (a heuristic
#'     method that tries to minimize the variance criterion under the
#'     restriction of equal group sizes). Anticlustering can be done via
#'     maximizinz the distance or the variance criterion.
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
#' anticlusters <- heuristic_anticlustering(features, preclusters)
#' # Check out results
#' plot(features, col = anticlusters, pch = 19)
#' tapply(features[, 1], anticlusters, mean)
#' tapply(features[, 1], anticlusters, sd)
#' tapply(features[, 2], anticlusters, mean)
#' tapply(features[, 2], anticlusters, sd)
#'
#' anticlusters <- heuristic_anticlustering(features, preclusters, objective = "variance")
#'
#' @references
#' H. Späth, “Anticlustering: Maximizing the variance criterion,”
#' Control and Cybernetics, vol. 15, no. 2, pp. 213–218, 1986.

heuristic_anticlustering <- function(features, clustering, objective = "distance",
                                     nrep = 100) {

  ## Some input handling
  features <- as.matrix(features) # if only one feature is passed
  if (!objective %in% c("distance", "variance"))
    stop("Argument objective must be 'distance' or 'variance'.")
  legal_number_of_clusters(features, clustering)

  ## Initialize variables
  n_elements <- nrow(features)
  n_preclusters <- length(unique(clustering))
  n_anticlusters <- n_elements / n_preclusters

  ## Store all data as a matrix for sorting. First column:
  ## Cluster affiliation; Second column: Original order of elements
  dat <- cbind(clustering, 1:n_elements, features)
  ## Order by precluster to ensure that each element from the same
  ## precluster is assigned to a different anticluster
  dat <- sort_by_col(dat, 1)

  ## Start optimizing
  best_obj <- -Inf

  for (i in 1:nrep) {
    anticlusters <- replicate_sample(n_preclusters, n_anticlusters)
    cur_obj <- get_objective(dat[, -(1:2)], anticlusters, objective)
    if (cur_obj > best_obj) {
      best_assignment <- anticlusters
      best_obj <- cur_obj
    }
  }
  ## Return anticluster assignment in original order
  dat[, 1] <- best_assignment
  dat <- sort_by_col(dat, 2)
  return(dat[, 1])
}

## Random anticluster assignment, replicated per precluster; called
## from within `heuristic_anticlustering`
replicate_sample <- function(times, N) {
  c(replicate(times, sample(N)))
}


#' Edit distances of neighbours
#'
#' Based on a preclustering, distances between items of the same cluster
#' are set to a large negative value. By considering a preclustering of
#' items, very similar items will not be assigned to the same group when
#' the fixed distance object is used to create the ILP formulation of
#' the item assignment instance.
#'
#' @param distances A distance object or matrix of
#'     between-item-distances.
#' @return A vector representing the group assignment; objects that are
#'     part of the same group as indicated by this vector are assigned a
#'     new distance.
#' @param value The value that is assigned to the fixed distances.
#'     Defaults to -1,000,000 currently.
#'
#' @return A distance object containing the fixed distances. For items
#'     that are part of the same cluster (as specified in the
#'     `data.frame` `assignment`), the between-item-distances are set to
#'     -1000000. This will have to be replaced by a theoretically-sound
#'     value; -1000000 is just a hack that will work in the present
#'     applications.
#'
#' @importFrom stats as.dist dist
#' @export

edit_distances <- function(distances, assignment, value = -1000000) {
  n_groups <- length(unique(assignment))
  distances <- as.matrix(distances)
  for (i in 1:n_groups) {
    items <- which(assignment == i) ## which items are in the i'th cluster
    ## two for-loops to tap all pairwise distances that may require
    ## fixing (a lot of unnecessary iterations, probably)
    for (j in 1:(length(items) - 1)) {
      for (t in 2:length(items)) {
        distances[items[j], items[t]] <- value
        distances[items[t], items[j]] <- value
      }
    }
  }
  return(as.dist(distances))
}
