
# Anticlustering based on a heuristic
#
# @param features A data.frame, matrix or vector representing the
#     features that are used.
# @param clustering A vector representing the preclustering of
#     elements. See Details.
# @param objective The objective to be maximized, either "distance"
#     (default) or "variance".
# @param nrep The number of repetitions tried when assigning elements
#     to anticlusters.
#
# @return A vector representing the anticlustering.
#
# @details The heuristic approach to anticlustering forbids elements
#     belonging to the same precluster to be assigned to the same
#     anticluster. The preclustering should be accomplished by one of
#     the clustering functions, `equal_sized_cluster_editing` (an exact
#     method that minimizes distance criterion under the restriction of
#     equal group sizes) or `equal_sized_kmeans` (a heuristic method
#     that tries to minimize the variance criterion under the
#     restriction of equal group sizes). Anticlustering can be done via
#     maximizinz the distance or the variance criterion.
#
# @examples
#
# features <- matrix(rnorm(1000, 100, 15), ncol = 2)
# n_anticlusters <- 4
# # Precluster cases
# n_preclusters <- nrow(features) / n_anticlusters
# preclusters <- equal_sized_kmeans(features, n_preclusters)
# # Use preclustering as resticting information in anticlustering
# anticlusters <- heuristic_anticlustering(features, preclusters)
# # Check out results
# plot(features, col = anticlusters, pch = 19)
# tapply(features[, 1], anticlusters, mean)
# tapply(features[, 1], anticlusters, sd)
# tapply(features[, 2], anticlusters, mean)
# tapply(features[, 2], anticlusters, sd)
#
# anticlusters <- heuristic_anticlustering(features, preclusters, objective = "variance")
#
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

