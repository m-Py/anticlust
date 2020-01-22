
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
#' par(mfrow = c(1, 1))
#'
#' @references
#'
#' Papenberg, M., & Klau, G. W. (2019, October 30). Using anticlustering
#' to partition a stimulus pool into equivalent parts.
#' https://doi.org/10.31234/osf.io/3razc
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
