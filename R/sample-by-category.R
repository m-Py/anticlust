
#' Random sampling employing a categorical constraint
#' 
#' This function can be used to obtain a stratified split of a data set.
#'
#' @param categories A matrix or vector of one or more categorical variables.
#' @param K The number of groups that are returned. 
#'
#' @return A vector representing the sample each element was assigned to.
#' 
#' @details 
#' 
#' This function can be used to obtain a stratified split of a data set. 
#' Using this function is like calling \code{\link{anticlustering}} with 
#' argument `categories`, but without optimizing a clustering objective. The
#' categories are just evenly split between samples. Apart from the restriction 
#' that categories are balanced between samples, the split is random.
#'
#' @export
#' 
#' @examples 
#' 
#' data(schaper2019)
#' categories <- schaper2019$room
#' groups <- categorical_sampling(categories, K = 6)
#' table(groups, categories)
#' 
#' # Unequal sized groups
#' groups <- categorical_sampling(categories, K = c(24, 24, 48))
#' table(groups, categories)
#' 

categorical_sampling <- function(categories, K) {
  categories <- merge_into_one_variable(categories)
  N <- length(categories)
  validate_input(K, "K", objmode = "numeric", 
                 must_be_integer = TRUE, not_na = TRUE)

  cats <- data.frame(
    categories = categories,
    order = 1:N
  )
  cats <- sort_by_col(cats, "categories")
  # Unequal group sizes were requested if sum(K) = N
  if (sum(K) == N) {
    samples <- replicate_sample(K)
    does_not_fit <- any(sort(table(samples)) != sort(K))
    while (does_not_fit) {
      samples <- replicate_sample(K)
      does_not_fit <- any(sort(table(samples)) != sort(K))
    }
  } else {
    validate_input(K, "K", objmode = "numeric", len = 1, 
                   greater_than = 1, must_be_integer = TRUE, not_na = TRUE)
    init <- rep_len(1:K, N)
    samples <- unlist(by(init, cats$categories, sample_))
  }
  cats$samples <- samples
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

# Weighted sampling for unequal-sized groups
sample_weighted <- function(tab) {
  n_categories <- length(tab)
  sample(1:n_categories, prob = tab / min(tab), size = sum(tab), replace = TRUE)
}

sample_stuff <- function(tab) {
  k <- length(tab)
  proportions <- tab / min(tab)
  has_to_be_in <- floor(proportions)
  sample <- rep(1:k, has_to_be_in)
  add <- which(as.logical(rbinom(k, 1, proportions - has_to_be_in)))
  sort(c(sample, add))
}

replicate_sample <- function(tab) {
  N <- sum(tab) 
  k <- length(tab)
  samples <- unlist(replicate(N / sum(floor(tab / min(tab))), sample_stuff(tab)))
  samples[1:sum(tab)]
}
