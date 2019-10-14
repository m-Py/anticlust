
#' Generate exchange partners for anticlustering
#'
#' This function is used to generate a vector that is used as the
#' \code{categories} argument in the function \code{\link{anticlustering}}.
#' This function is meant to reduce the total number of exchange partners
#' (the default is *all* elements) so that the exchange method runs
#' faster.
#'
#' @param N The number of elements. Does not need to be passed if either
#' \code{categories}, \code{features} or \code{distance} is passed.
#' @param p The number of exchange partners per element.
#' @param categories A vector describing categories of elements.
#' @param similar Logical. Should similar items serve as exchange partners?
#' Defaults to \code{FALSE}.
#' @param features Optional data argument. If the exchange partners
#' are generated via a cluster analysis (that is if argument
#' \code{similar = TRUE}), data is needed describing the elements.
#'
#' @details
#' This function will be particularly interesting when generating exchange
#' partners within categories. For imbalanced data (that is, within each
#' category, the elements cannot be split into equal-sized parts of
#' size \code{p}, the function may return some elements with \code{p + 1}
#' exchange partners.
#'
#' In general, the function will try to balance out everything as well as
#' possible (this is also true for clustering and especially clustering
#' within categories), but it cannot always succeed. That is, sometimes
#' the function will throw an error if the argument \code{p} does not
#' fit with the arguments \code{similar = TRUE} (which request a cluster
#' analysis with clusters of equal size) or \code{categories}.
#'
#' @export
#'

generate_exchange_partners <- function(
  N = NULL,
  p,
  categories = NULL,
  similar = FALSE,
  features = NULL
) {
  if (argument_exists(N)) {
    return(sample(rep_len(1:p, N)))
  }
  if (argument_exists(categories) && similar == FALSE) {
    # some sorting etc
    categories <- merge_into_one_variable(categories)
    categories <- data.frame(1:length(categories), categories)
    categories <- sort_by_col(categories, 2)
    partners <- lapply(unique(categories[, 2]), function(x) {
      n <- sum(categories[, 2] == x)
      k <- p + 1
      sample(rep_len(1:(n/k), n))
    })
    categories$partners <- to_numeric(paste0(unlist(partners), categories[, 2]))
    return(sort_by_col(categories, 1)$partners)
  }
  if (!argument_exists(categories) && similar == TRUE) {
    K <- nrow(features) / (p + 1)
    print(K)
    return(balanced_clustering(features, K = K))
  }
  if (argument_exists(categories) && similar == TRUE) {
    K <- p + 1
    if (!preclustering_possible(categories, K)) {
      stop("Given your categories, the data cannot be split into clusters of size p + 1")
    }
    return(to_numeric(precluster_per_category(features, categories, K)))
  }
  stop("Pass one of the arguments N, categories, features to get a result")
}
