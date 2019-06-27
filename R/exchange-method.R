

#' Solve anticlustering using the modified exchange method
#'
#' @param data the data -- a N x N dissimilarity matrix or a N x M
#'     table of item features
#' @param K The number of anticlusters
#' @param obj_function the objective function (to compute ACE or
#'     K-Means criterion). Takes as first argument a cluster assignment
#'     and as second argument the data set `data`.#'
#' @param categories A vector representing categorical constraints
#'
#' @return The anticluster assignment
#'
#' @noRd
#'
#' @export
#'


exchange_method <- function(data, K, obj_function, categories = NULL) {
  ## generate a legal cluster assignment, satisfying the categorical constraints
  if (!is.null(categories)) {
    clusters <- heuristic_anticlustering(data, K, NULL, "distance",
                                         1, NULL, categories, FALSE, NULL, ncores = NULL)
  } else {
    clusters <- rep_len(1:K, length.out = nrow(data))
  }

  N <- nrow(data)
  best_total <- obj_function(clusters, data)
  for (i in 1:N) {
    # cluster of current item
    tmp_group <- clusters[i]
    # are there categorical variables?
    if (!is.null(categories)) {
      # only exchange within the same group
      allowed_partner <- categories == categories[i]
    } else {
      allowed_partner <- rep(TRUE, nrow(data)) # no constraint
    }
    exchange_partners <- (clusters != tmp_group) & allowed_partner
    # items in other clusters
    exchange_partners <- (1:N)[exchange_partners]
    # container to store objectives associated with each exchange of item i:
    comparison_objectives <- rep(NA, N)
    for (j in exchange_partners) {
      ## Swap all pairs of items and check out objective
      tmp_clusters <- clusters
      tmp_clusters[i] <- tmp_clusters[j]
      tmp_clusters[j] <- tmp_group
      comparison_objectives[j] <- variance_objective_(tmp_clusters, data)
    }
    ## Do the swap if an improvement occured
    best_this_round <- max(comparison_objectives, na.rm = TRUE)
    if (best_this_round > best_total) {
      # which element has to be swapped
      swap <- which(comparison_objectives == best_this_round)[1]
      # swap the elements
      clusters[i] <- clusters[swap]
      clusters[swap] <- tmp_group
      # update best solution
      best_total <- best_this_round
    }
  }
  clusters
}
