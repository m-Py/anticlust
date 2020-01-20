
#' Random sampling employing a categorical constraint
#'
#' @param categories A vector of categories
#' @param K The number of anticlusters
#'
#' @return A random shuffling of the anticlusters that balances out
#'     the categories between samples
#'
#' @noRd

categorical_sampling <- function(categories, K) {
  N <- length(categories)
  cats <- data.frame(
    categories = categories,
    order = 1:N
  )
  cats <- sort_by_col(cats, "categories")
  cats$samples <- unlist(by(rep_len(1:K, N), cats$categories, sample_))
  sort_by_col(cats, "order")$samples
}

# sample function must be redefined for categorical sampling just in case
# there is a category with only 1 member
sample_ <- function(x, ...) {
  if (length(x) == 1) {
    return(x)
  }
  sample(x, ...)
}
