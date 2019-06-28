
#' Wrapper for the exchange method algorithm below
#'
#' This function takes the input and determines the objective function;
#' i.e., just some preprocessing for the algorithm. (
#'
#' For the arguments, see \code{anticlustering}
#'
#' @noRd
#'

exchange_method <- function(features, distances, K, objective,
                            categories) {
  categories <- merge_into_one_variable(categories) # may be NULL
  obj_function <- get_objective_function(features, distances, objective)
  if (argument_exists(features)) {
    data <- features
  } else {
    data <- distances
  }
  exchange_method_(data, K, obj_function, categories)
}

#' Solve anticlustering using the modified exchange method
#'
#' @param data the data -- a N x N dissimilarity matrix or a N x M
#'     table of item features
#' @param K The number of anticlusters
#' @param obj_function the objective function (to compute ACE or
#'     K-Means criterion). Takes as first argument a cluster assignment
#'     and as second argument the data set `data`.
#' @param categories A vector representing categorical constraints
#'
#' @return The anticluster assignment
#'
#' @noRd
#'
#'

exchange_method_ <- function(data, K, obj_function, categories = NULL) {
  if (!is.null(categories)) {
    # An initial cluster assignment satisfying the categorical constraints
    clusters <- random_sampling(data, K, NULL, "distance",
                                1, NULL, categories, FALSE,
                                NULL, ncores = NULL)
  } else {
    clusters <- rep_len(1:K, length.out = nrow(data))
  }

  N <- nrow(data)
  best_total <- obj_function(clusters, data)
  for (i in 1:N) {
    # cluster of current item
    group_i <- clusters[i]
    # are there categorical variables?
    if (!is.null(categories)) {
      # only exchange within the same group
      allowed_partner <- categories == categories[i]
    } else {
      allowed_partner <- rep(TRUE, nrow(data)) # no constraint
    }
    exchange_partners <- (clusters != group_i) & allowed_partner
    # items in other clusters
    exchange_partners <- (1:N)[exchange_partners]
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, N)
    for (j in exchange_partners) {
      ## Swap item i with all legal exchange partners and check out objective
      tmp_clusters <- clusters
      tmp_clusters[i] <- tmp_clusters[j]
      tmp_clusters[j] <- group_i
      comparison_objectives[j] <- variance_objective_(tmp_clusters, data)
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives, na.rm = TRUE)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- which(comparison_objectives == best_this_round)[1]
      # swap the elements
      clusters[i] <- clusters[swap]
      clusters[swap] <- group_i
      # update best solution
      best_total <- best_this_round
    }
  }
  clusters
}
