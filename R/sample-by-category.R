
#' Random sampling employing a categorical constraint
#' 
#' This function can be used to obtain a stratified split of a data set.
#'
#' @param categories A vector of a categorical variable.
#' @param K The number of groups that are returned. 
#'
#' @return A vector representing the sample each element was assigned to.
#' 
#' @details 
#' 
#' This function can be used to obtain a stratified split of a data set. 
#' Using this function is like calling a\code{\link{anticlustering}}` with 
#' argument `categories` where no optimization is conducted; the categories are 
#' just evenly split between samples. Apart from the restriction that categories 
#' are balanced between samples, the split is random.
#'
#' @export

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
