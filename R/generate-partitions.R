
#' Generate all partitions of same cardinality
#'
#' @param K How many partitions
#' @param N The total N. \code{K} has to be dividble
#'     by distinct.
#' @param generate_permutations If TRUE, all permutations are returned,
#'     resulting in duplicate partitions.
#'
#' @param A list of all partitions.
#'
#' @examples
#'
#' start <- Sys.time()
#' partitions <- generate_partitions(2, 10)
#' intermediate <- Sys.time()
#' permutations <- generate_partitions(2, 10, TRUE)
#' end <- Sys.time()
#' intermediate - start
#' end - intermediate
#'
#' @noRd
#'

generate_partitions <- function(K, N, generate_permutations = FALSE) {
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
