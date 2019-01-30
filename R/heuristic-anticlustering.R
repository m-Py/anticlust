
# Anticlustering based on a heuristic
#
# @param features A data.frame, matrix or vector representing the
#     features that are used.
# @param clustering A vector representing the preclustering of
#     elements. See Details.
# @param objective The objective to be maximized, either "distance"
#     (default) or "variance".
# @param method The method used to find the best objective. Can be "rnd"
#     for repeated random sampling or "sa" for simulated annealing.
# @param nrep The number of repetitions tried when assigning elements
#     to anticlusters when the method is "rnd"
#
# @return A vector representing the anticlustering.
#
# @details
#     The heuristic approaches to anticlustering forbids elements
#     belonging to the same precluster to be assigned to the same
#     anticluster. The preclustering should be accomplished by one of
#     the clustering functions, `equal_sized_cluster_editing` (an exact
#     method that minimizes distance criterion under the restriction of
#     equal group sizes) or `equal_sized_kmeans` (a heuristic method
#     that tries to minimize the variance criterion under the
#     restriction of equal group sizes). Anticlustering can be done via
#     maximizinz the distance or the variance criterion.
#     The simulated annealing approach also considers the preclusters
#     as the neighborhood in the exchange candidate generation. Only
#     exchanges between elements of the same precluster are allowed
#     in the simulated annealing approach.
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
                                     method = "rnd", nrep = 100) {

  ## Some input handling
  if (!objective %in% c("distance", "variance"))
    stop("Argument objective must be 'distance' or 'variance'.")
  legal_number_of_clusters(features, clustering)

  ## Store all data as a matrix for sorting. First column:
  ## Cluster affiliation; Second column: Original order of elements
  dat <- cbind(clustering, 1:nrow(features), features)
  ## Order by precluster to ensure that each element from the same
  ## precluster is assigned to a different anticluster
  dat <- sort_by_col(dat, 1)

  if (method == "rnd")
    best_assignment <- random_sampling(dat, clustering, objective, nrep)
  else if (method == "sa")
    best_assignment <- simulated_annealing(dat, clustering, objective)
  else
    stop("Method must be 'rnd' or 'sa'")

  ## Return anticluster assignment in original order
  dat[, 1] <- best_assignment
  dat <- sort_by_col(dat, 2)
  return(dat[, 1])
}

## Simulated annealing approach to finding the best objective considering
## preclustering restrictions

simulated_annealing <- function(dat, clustering, objective) {
  ## Initialize variables
  n_elements <- nrow(dat)
  n_preclusters <- length(unique(clustering))
  n_anticlusters <- n_elements / n_preclusters
  ## Initial parameter values
  init <- replicate_sample(n_preclusters, n_anticlusters)

  ## Wrap objective function so it takes only one parameter
  objective_fun <- function(clusters) {
    return(get_objective(dat[, -c(1, 2)], clusters, objective))
  }

  ## Use `optim` for simulated annealing
  return(optim(init, objective_fun, next_candidate, method = "SANN",
               control = list(maxit = 10000, temp = 2000, REPORT = 500,
                              fnscale = -1, tmax = 20))$par)
}

## Random sampling approach to finding the best objective considering
## preclustering restrictions
random_sampling <- function(dat, clustering, objective, nrep) {
  ## Initialize variables
  n_elements <- nrow(dat)
  n_preclusters <- length(unique(clustering))
  n_anticlusters <- n_elements / n_preclusters
  ## Start optimizing
  best_obj <- -Inf
  for (i in seq_along(nrep)) {
    anticlusters <- replicate_sample(n_preclusters, n_anticlusters)
    cur_obj <- get_objective(dat[, -(1:2)], anticlusters, objective)
    if (cur_obj > best_obj) {
      best_assignment <- anticlusters
      best_obj <- cur_obj
    }
  }
  return(best_assignment)
}

## Random anticluster assignment, replicated per precluster; called
## from within `heuristic_anticlustering`
replicate_sample <- function(times, N) {
  c(replicate(times, sample(N)))
}


## Generate a candidate move for simulated annealing anticlustering.
## This function only works as expected when the anticlusters are sorted
## by precluster, as is the case when it is called from within
next_candidate <- function(anticlusters) {
  n_anticlusters <- length(unique(anticlusters))
  n_preclusters <- length(anticlusters) / n_anticlusters
  preclusters <- rep(1:n_preclusters, each = n_anticlusters)
  ## select a random precluster
  rndclus <- sample(n_preclusters, 1)
  ## within this precluster, which elements should be swapped?
  changepoints <- sample(n_anticlusters, size = 2, replace = FALSE)
  ## Some ugly code for swapping the two values
  tmp <- anticlusters[preclusters == rndclus][changepoints[1]]
  anticlusters[preclusters == rndclus][changepoints[1]] <- anticlusters[preclusters == rndclus][changepoints[2]]
  anticlusters[preclusters == rndclus][changepoints[2]] <- tmp
  return(anticlusters)
}
