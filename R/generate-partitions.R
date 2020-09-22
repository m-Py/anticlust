
#' Generate all partitions of same cardinality
#'
#' @param N The total N. \code{K} has to be dividble
#'     by \code{N}.
#' @param K How many partitions
#' @param generate_permutations If TRUE, all permutations are returned,
#'     resulting in duplicate partitions.
#'
#' @return A list of all partitions (or permutations if
#' \code{generate_permutations} is \code{TRUE}).
#'
#' @details
#'
#' In principle, anticlustering can be solved to optimality by
#' generating all possible partitions of N items into K groups.
#' The example code below illustrates how to do this.
#' However, this approach only works for small N because the
#' number of partitions grows exponentially with N.
#'
#' The partition c(1, 2, 2, 1)
#' is the same as the partition c(2, 1, 1, 2) but they correspond
#' to different permutations of the elements [1, 1, 2, 2]. If the argument
#' \code{generate_permutations} is \code{TRUE}, all permutations are
#' returned. To solve balanced anticlustering exactly, it is sufficient
#' to inspect all partitions while ignoring duplicated permutations.
#'
#' @examples
#'
#' ## Generate all partitions to solve k-means anticlustering
#' ## to optimality.
#'
#' N <- 14
#' K <- 2
#' features <- matrix(sample(N * 2, replace = TRUE), ncol = 2)
#' partitions <- generate_partitions(N, K)
#' length(partitions) # number of possible partitions
#'
#' ## Create an objective function that takes the partition
#' ## as first argument (then, we can use sapply to compute
#' ## the objective for each partition)
#' var_obj <- function(clusters, features) {
#'   variance_objective(features, clusters)
#' }
#'
#' all_objectives <- sapply(
#'   partitions,
#'   FUN = var_obj,
#'   features = features
#' )
#'
#' ## Check out distribution of the objective over all partitions:
#' hist(all_objectives) # many large, few low objectives
#' ## Get best k-means anticlustering objective:
#' best_obj <- max(all_objectives)
#' ## It is possible that there are multiple best solutions:
#' sum(all_objectives == best_obj)
#' ## Select one best partition:
#' best_anticlustering <- partitions[all_objectives == best_obj][[1]]
#' ## Look at mean for each partition:
#' by(features, best_anticlustering, function(x) round(colMeans(x), 2))
#'
#'
#' ## Get best k-means clustering objective:
#' min_obj <- min(all_objectives)
#' sum(all_objectives == min_obj)
#' ## Select one best partition:
#' best_clustering <- partitions[all_objectives == min_obj][[1]]
#'
#' ## Plot minimum and maximum objectives:
#' user_par <- par("mfrow")
#' par(mfrow = c(1, 2))
#' plot_clusters(
#'   features,
#'   best_anticlustering,
#'   illustrate_variance = TRUE,
#'   main = "Maximum variance"
#' )
#' plot_clusters(
#'   features,
#'   best_clustering,
#'   illustrate_variance = TRUE,
#'   main = "Minimum variance"
#' )
#' par(mfrow = user_par)
#'
#' @references
#'
#' Papenberg, M., & Klau, G. W. (2020). Using anticlustering to partition 
#' data sets into equivalent parts. Psychological Methods. Advance Online 
#' Publication. https://doi.org/10.1037/met0000301.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'

generate_partitions <- function(N, K, generate_permutations = FALSE) {
  if (N %% K != 0) {
    stop("K must be a divider of N.")
  }
  partitions <- 1:K
  init <- rep(partitions, each = N / K)
  anticlusters <- init
  all_ <- list(anticlusters)
  i <- 2
  while (!all(next_permutation(anticlusters) == init)) {
    anticlusters <- next_permutation(anticlusters)
    if (anticlusters[1] == max(anticlusters) &&
        generate_permutations == FALSE) {
      break
    }
    if (!is.unsorted(match(partitions, anticlusters)) ||
        generate_permutations == TRUE) {
      all_[[i]] <- anticlusters
      i <- i + 1
    }
  }
  return(all_)
}

#' Get the next permutation 
#'
#' Each permutation is computed on basis of a passed permutation where
#' lexicographic ordering of the permutations is used to determine the
#' next "higher" permutation. This is an adaption of the
#' `next_permutation` function in C++.
#'
#' @param permutation A vector of elements.
#'
#' @return The next higher permutation of the elements in vector
#'     `permutation` with regard to its lexicographic ordering.
#'
#' 
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @noRd
#' 
next_permutation <- function(permutation) {
  n    <- length(permutation)
  last <- permutation[n]
  i    <- n
  while(last <= permutation[i-1]) {
    last <- permutation[i-1]
    i    <- i - 1
    ## if lexicographic order is already at the maximum:
    if (i-1 == 0) return(sort(permutation))
  }
  ## this algorithm divides the input in a head and a tail; the tail
  ## is monotonically decreasing
  head <- permutation[1:(i-1)]
  tail <- permutation[i:length(permutation)]
  ## which element in the tail is the smallest element that is larger
  ## than the last element in the head?
  larger_values  <- tail[tail > head[length(head)]]
  ## last element of the head:
  final_head <- head[length(head)]
  ## replace last element of head by smallest larger value in tail
  head[length(head)] <- min(larger_values)  
  ## replace smallest larger value in tail by final head element
  tail[max(which(tail == min(larger_values)))] <- final_head
  ## reverse tail before returning
  return(c(head, rev(tail)))
}
