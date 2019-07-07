
#' Wrapper for the exchange method algorithm below
#'
#' This function takes the input and determines the objective function;
#' i.e., just some preprocessing for the algorithm.
#'
#' For the arguments, see \code{anticlustering}
#'
#' @noRd
#'

exchange_method <- function(features, distances, K, obj_function,
                            categories, preclusters) {
  if (argument_exists(features)) {
    data <- features
  } else {
    data <- distances
  }

  ## Generate an initial clustering from which the
  ## exchange is conducted. Default case is random initiation:
  clusters <- sample(rep_len(1:K, length.out = nrow(data)))

  ## The problem is: An initial assingnment is needed that potentially
  ## satisfies constraints (preclustering and/or categorical).
  if (!is.null(preclusters)) {
    clusters <- random_sampling(features, K, preclusters, obj_function,
                                1, distances, NULL, FALSE,
                                NULL, NULL)
  }
  if (!is.null(categories)) {
    ## Categorical constraints have higher priority and overwrite
    ## preclustering constraints. However, the exchange algorithm
    ## below still tries to adhere to the preclustering constraints as
    ## well as possible
    clusters <- random_sampling(data, K, NULL, obj_function,
                                1, distances, categories, FALSE,
                                NULL, NULL)
  }

  exchange_method_(data, clusters, obj_function, categories, preclusters)
}

#' Solve anticlustering using the modified exchange method
#'
#' @param data the data -- a N x N dissimilarity matrix or a N x M
#'     table of item features
#' @param clusters An initial cluster assignment
#' @param obj_function the objective function (to compute ACE or
#'     K-Means criterion). Takes as first argument a cluster assignment
#'     and as second argument the data set `data`.
#' @param categories A vector representing categorical constraints
#' @param preclusters A vector representing preclustering constraints
#'
#' @return The anticluster assignment
#'
#' @noRd
#'
#'

exchange_method_ <- function(data, clusters, obj_function, categories, preclusters) {
  N <- nrow(data)
  best_total <- obj_function(clusters, data)
  for (i in 1:N) {
    # cluster of current item
    group_i <- clusters[i]
    # are there categorical variables?
    if (!is.null(categories)) {
      # only exchange within the same group
      allowed_category <- categories == categories[i]
    } else {
      allowed_category <- rep(TRUE, nrow(data)) # no constraint
    }
    if (!is.null(preclusters)) {
      allowed_precluster <- preclusters == preclusters[i]
    } else {
      allowed_precluster <- rep(TRUE, nrow(data)) # no constraint
    }

    ## Feasible exchange partners:
    # (a) are in different anticluster
    # (b) are in the same precluster
    # (c) have the same category
    exchange_partners <- (clusters != group_i) & allowed_category & allowed_precluster
    ## Ensure there are more than zero exchange partners:
    if (sum(exchange_partners) == 0) {
      exchange_partners <- (clusters != group_i) & allowed_category
    }
    # items in other clusters
    exchange_partners <- (1:N)[exchange_partners]
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, length(exchange_partners))
    for (j in seq_along(exchange_partners)) {
      ## Swap item i with all legal exchange partners and check out objective
      tmp_clusters <- clusters
      tmp_clusters[i] <- tmp_clusters[exchange_partners[j]]
      tmp_clusters[exchange_partners[j]] <- group_i
      comparison_objectives[j] <- variance_objective_(tmp_clusters, data)
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- exchange_partners[comparison_objectives == best_this_round][1]
      # swap the elements
      clusters[i] <- clusters[swap]
      clusters[swap] <- group_i
      # update best solution
      best_total <- best_this_round
    }
  }
  clusters
}
