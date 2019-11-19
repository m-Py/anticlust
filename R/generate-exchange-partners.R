
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
#'
#' This function is used to define which elements serve as exchange
#' partners when calling the anticlustering exchange method. Its output
#' is passed to the function \code{\link{anticlustering}} as the argument
#' \code{categories}.
#'
#' This function is usually used to reduce the number of exchange
#' partners to speed up the exchange method. By default, the exchange
#' method swaps each element with each other element (unless the other
#' element is already part of the same cluster). When passing a
#' categorical vector, only the member of the same category serve as
#' exchange partner. This function is particularly useful to ensure
#' that not even all elements from the same category are used as exchange
#' partners, but only a subset.
#'
#' Note that this function may sometimes run into problems when dealing
#' with imbalanced data (that is, within each category, the elements
#' cannot be split into equal-sized parts of size \code{p + 1}).
#' Therefore, the function may return more than \code{p} exchange
#' partners for some elements. When \code{similar = TRUE}, such impossible
#' imbalances will lead to an error message.
#'
#' If the argument \code{similar} is \code{TRUE}, the exchange partners
#' are chosen in such a way that similar items are exchange partners
#' for each other (following the general preclustering logic, see
#' \code{\link{anticlustering}}) Similar items are grouped using a call
#' to the function \code{\link{balanced_clustering}}. When combining the
#' arguments \code{similar = TRUE} and \code{categories}, within each
#' category a cluster analysis is conducted to group similar items.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' generate_exchange_partners(N = 6, p = 1)
#' # impossible requirement: each of 19 elements has exactly 3 partners:
#' table(generate_exchange_partners(N = 19, p = 3))
#'
#' ## On the schaper2019 data set
#' data(schaper2019)
#' # Enforce 3 exchange partners per element - faster than using each item
#' start <- Sys.time()
#' anticlustering(
#'   schaper2019[, 3:6],
#'   K = 3,
#'   categories = generate_exchange_partners(categories = schaper2019$room, p = 3)
#' )
#' Sys.time() - start
#'
#' # Without these restrictions:
#' start <- Sys.time()
#' anticlustering(
#'   schaper2019[, 3:6],
#'   K = 3,
#'   categories = schaper2019$room
#' )
#' Sys.time() - start
#'
#' ## Pass a data set and p:
#' table(generate_exchange_partners(features = rnorm(999), p = 2))
#'

generate_exchange_partners <- function(
  N = NULL,
  p,
  categories = NULL,
  similar = FALSE,
  features = NULL
) {
  if (argument_exists(N)) {
    return(sample_partners(N, p))
  }
  if (argument_exists(categories) && similar == FALSE) {
    # some sorting etc
    categories <- merge_into_one_variable(categories)
    categories <- data.frame(1:length(categories), categories)
    categories <- sort_by_col(categories, 2)
    partners <- lapply(unique(categories[, 2]), function(x) {
      n <- sum(categories[, 2] == x)
      sample_partners(n, p)
    })
    categories$partners <- to_numeric(paste0(unlist(partners), categories[, 2]))
    return(sort_by_col(categories, 1)$partners)
  }
  if (!argument_exists(categories) && similar == TRUE) {
    K <- p + 1
    return(imbalanced_preclustering(features, K = K))
  }
  if (argument_exists(categories) && similar == TRUE) {
    K <- p + 1
    return(to_numeric(precluster_per_category(features, categories, K)))
  }
  if (argument_exists(features)) {
    return(sample_partners(NROW(features), p))
  }
  stop("Pass one of the arguments N, categories, features to get a result")
}

sample_partners <- function(N, p) {
  k <- p + 1
  sample(rep_len(1:(N/k), N))
}
