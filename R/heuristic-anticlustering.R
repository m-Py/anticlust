
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
# @param preclustering Boolean, should a preclustering be conducted
#     before anticlusters are created. Defaults to TRUE and it is
#     advised to keep it that way for the heuristic methods.
#     `preclustering` = FALSE is mainly implemented to test against the
#     option `preclustering` = TRUE
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
                                     method = "rnd", nrep, preclustering = TRUE) {

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
    best_assignment <- random_sampling(dat, clustering, objective, nrep, preclustering)
  else if (method == "sa")
    best_assignment <- simulated_annealing(dat, clustering, objective, preclustering)
  else
    stop("Method must be 'rnd' or 'sa'")

  ## Return anticluster assignment in original order
  dat[, 1] <- best_assignment
  dat <- sort_by_col(dat, 2)
  return(dat[, 1])
}

## Simulated annealing approach to finding the best objective considering
## preclustering restrictions. Relevant for this method is the
## next_candidate function below
simulated_annealing <- function(dat, clustering, objective, preclustering = TRUE) {
  ## Initialize variables
  n_elements <- nrow(dat)
  n_preclusters <- length(unique(clustering))
  n_anticlusters <- n_elements / n_preclusters
  ## Initial parameter values; the preclustering ensures a good initial
  ## state for the simulated annealing approach. Uses a random initial
  ## state when preclustering is not considered.
  if (preclustering == FALSE) {
    init <- rep(1:n_anticlusters, n_preclusters)
    init <- sample(anticlusters)
  } else if (preclustering == TRUE) {
    init <- replicate_sample(n_preclusters, n_anticlusters)
  }

  ## Wrap objective function so it only takes the anticluster
  ## affiliation as parameter (data is included from outside the function)
  objective_fun <- function(clusters) {
    return(get_objective(dat[, -c(1, 2)], clusters, objective))
  }

  ## Use `optim` for simulated annealing
  return(optim(init, objective_fun, next_candidate, method = "SANN",
               preclustering = preclustering,
               control = list(maxit = 10000, temp = 2000, REPORT = 500,
                              fnscale = -1, tmax = 20))$par)
}


## Generate a candidate move for simulated annealing anticlustering.
## This function only works as expected when the anticlusters are sorted
## by precluster, as is the case when it is called from within the
## `simulated_annealing` function above.
next_candidate <- function(anticlusters, preclustering = TRUE) {
  n_anticlusters <- length(unique(anticlusters))
  n_preclusters <- length(anticlusters) / n_anticlusters
  preclusters <- rep(1:n_preclusters, each = n_anticlusters)

  ## 1. Candidate generation without preclustering. Change the
  ## anticlusters of two entirely random elements
  if (preclustering == FALSE) {
    changepoints <- sample(anticlusters, size = 2, replace = FALSE)
    tmp <- anticlusters[changepoints[1]]
    anticlusters[changepoints[1]] <- anticlusters[changepoints[2]]
    anticlusters[changepoints[2]] <- tmp
    return(anticlusters)
  }

  ## 2. Next candidate generation while respecting preclustering

  ## Select a random precluster
  rndclus <- sample(n_preclusters, 1)
  ## Within this precluster, which elements should be swapped into
  ## the respective other anticluster?
  changepoints <- sample(n_anticlusters, size = 2, replace = FALSE)
  ## Some ugly code for swapping anticluster affiliation of the two elements
  tmp <- anticlusters[preclusters == rndclus][changepoints[1]]
  anticlusters[preclusters == rndclus][changepoints[1]] <- anticlusters[preclusters == rndclus][changepoints[2]]
  anticlusters[preclusters == rndclus][changepoints[2]] <- tmp
  return(anticlusters)
}


## Random sampling approach to finding the best objective considering
## preclustering restrictions
random_sampling <- function(dat, clustering, objective, nrep, preclustering = TRUE) {
  ## Initialize variables
  n_elements <- nrow(dat)
  n_preclusters <- length(unique(clustering))
  n_anticlusters <- n_elements / n_preclusters
  ## Start optimizing
  best_obj <- -Inf
  for (i in seq_along(nrep)) {
    ## 1. Random sampling without preclustering restrictions
    if (preclustering == FALSE) {
      anticlusters <- rep(1:n_anticlusters, n_preclusters)
      anticlusters <- sample(anticlusters)
    ## 2. Include preclustering restrictions
    } else if (preclustering == TRUE) {
      anticlusters <- replicate_sample(n_preclusters, n_anticlusters)
    }
    cur_obj <- get_objective(dat[, -(1:2)], anticlusters, objective)
    if (cur_obj > best_obj) {
      best_assignment <- anticlusters
      best_obj <- cur_obj
    }
  }
  return(best_assignment)
}

## Random anticluster assignment, replicated per precluster; called
## from within `random_sampling` and `simulated_annealing`
replicate_sample <- function(times, N) {
  c(replicate(times, sample(N)))
}

