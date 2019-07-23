
#' Generate an initial cluster assignment for anticlustering
#'
#' This function returns a vector that can be used as the \code{K} argument
#' in the function \code{\link{anticlustering}}. Some entries in this
#' vector may be NA to enable subset selection procedures (i.e., not each
#' item is assigned to a set).
#'
#' @param N The total number of items
#' @param K The number of subsets to be generated
#' @param n Either: a number indicating how many elements are selected
#'     per subset. Or: A vector of length K. In the latter case, entries
#'     represents the n per subset.
#' @param groups A vector of length N, alternative argument to N and K.
#'     Represents the group for an element.
#'
#' @return The initialized cluster for each element
#'
#' @details
#'
#' If the sum of `n` per group is lower than `N`, the returned vector
#' will include NAs. If this vector is used as the \code{K} argument
#' for the \code{\link{anticlustering}}, \code{anticlustering} will
#' also output a vector that contains NAs (that is a permutation of the
#' input vector \code{K}). This way, it is possible to select only a
#' subset of the input items. It is also possible to create sets of
#' different sizes.
#'
#' @examples
#'
#' # Example of how to use the function `anticlustering` to create two
#' # subset of stimuli from the schaper2019 data set. Two sets of n = 20
#' # each are created, one with kitchen items only and one with bathroom
#' # items only. The sets are parallelized on the means and standard
#' # deviations of four features.
#'
#' head(schaper2019)
#' features <- schaper2019[, 3:6]
#'
#' K <- initialize_K(groups = schaper2019$room, n = c(20, 20))
#' groups <- anticlustering(
#'   features,
#'   K = K,
#'   categories = schaper2019$room,
#'   objective = mean_sd_obj # this objective function makes most sense for subset selection
#' )
#' # Note that this subset selection based on two groups is
#' # no longer an anticlustering method (that would divide one pool of
#' # items into subsets)
#'
#' table(K, schaper2019$room)
#' # Compare feature means by room
#' by(features, groups, function(x) round(colMeans(x), 2))
#' # Compare standard deviations by room
#' by(features, groups, function(x) round(apply(x, 2, sd), 2))
#'
#'
#' # Create sets of different size:
#' K <- initialize_K(n = c(48, 24, 24))
#' anticlusters <- anticlustering(
#'   features,
#'   K = K,
#'   objective = "distance"
#' )
#'
#' # Compare feature means by room
#' table(anticlusters)
#'
#'
#' @export

initialize_K <- function(N = NULL, K = NULL, n, groups = NULL) {
  if (argument_exists(N) && argument_exists(groups)) {
    stop("do not pass both `N` and `groups` argument")
  }
  if (argument_exists(K) && argument_exists(groups)) {
    stop("do not pass both `K` and `groups` argument")
  }
  if (!argument_exists(K) && !argument_exists(groups)) {
    if (!argument_exists(N)) {
      if (length(n) == 1) {
        stop("If only the argument `n` is passed, `n` must be of length > 1")
      }
      return(sample(rep(1:length(n), n)))
    } else {
      return(initialize_K(N = N, K = length(n), n))
    }
  }

  if (argument_exists(N) && argument_exists(K)) {
    groups <- rep_len(1:K, N)
    N <- rep(N / K, K) # N per group
  }
  else if (argument_exists(groups)) {
    groups <- as.numeric(as.factor(groups))
    N <- table(groups)
    K <- length(N)
  }
  if (length(n) == 1) {
    n <- rep(n, K)
  }

  slots <- rep(99, sum(N))
  for (i in 1:K) {
    ## how many NAs must be filled:
    fill_NA <-  N[i] - n[i]
    real_slots <- N[i] - fill_NA
    slots[groups == i] <- sample(c(rep(i, real_slots), rep(NA, fill_NA)))
  }
  slots
}
