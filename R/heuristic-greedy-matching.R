
#' Greedy matching
#'
#' Finds pairs of elements that are similar. Can be used as a
#' preclustering step for heuristic anticluster editing when K = 2.
#'
#' @param distances A N x N matrix representing the
#'     pairwise dissimilarities between N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}).
#'
#' @details
#'
#' This function realizes a greedy matching approach given pairwise
#' dissimilarities of N elements. The algorithm works as follows:
#'
#' 1. Determine the two closest elements with regard to the
#'    pairwise dissimilarity
#'
#' 2. Assign these two elements to the same cluster
#'
#' 3. Remove these two elements from the pool
#'
#' 4. Repeat 1-3 until no elements are remaining
#'
#' @return A vector representing which elements were matched.
#'
#' @export
#'
#' @examples
#'
#' # Compare exact and greedy matching
#' N <- 10
#' n_features <- 2
#' features <- matrix(rnorm(N * n_features), ncol = n_features)
#' c1 <- greedy_matching(dist(features))
#' c2 <- balanced_clustering(features, K = N / 2, method = "exact")
#' par(mfrow = c(1, 2))
#' plot_clusters(features, c1, within_connection = TRUE, cex = 2.5)
#' plot_clusters(features, c2, within_connection = TRUE, cex = 2.5)
#'

greedy_matching <- function(distances) {

  ## Prepare the distance matrix
  distances <- as.matrix(distances)
  diag(distances) <- Inf
  N <- nrow(distances)
  if (N %% 2 != 0) {
    stop("The number of elements has to be a divider of 2.")
  }

  ## Repeatedly select the closest items and assign them to a cluster
  anticlusters <- rep(NA, N)
  for (i in 1:(N/2)) {
    ## Select the closest pair (if the smallest distance occurs more
    ## than once, select the first pair)
    pair <- which(distances == min(distances), arr.ind = TRUE)[1, ]
    ## Randomly assign the two stimuli to two different clusters:
    anticlusters[pair] <- i
    ## Remove used stimuli from the pool
    distances[pair, ] <- Inf
    distances[, pair] <- Inf
  }
  ## Assert that the output has the correct structure:
  if (!(all(table(anticlusters) == 2))) {
    stop("Something went wrong: not the same set sizes")
  }
  if (any(is.na(anticlusters))) {
    stop("Something went wrong: NA in output")
  }
  return(anticlusters)
}
