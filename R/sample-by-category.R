
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
#' @importFrom stats rbinom
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
#' # Heavily unequal sized groups, is harder to balance the groups
#' groups <- categorical_sampling(categories, K = c(51, 19, 26))
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
  if (sum(K) == N && length(K) > 1) { # Unequal group sizes were requested
    gdc <- gcd_set(K)
    if (gdc > 1) {
      # if the GDC of all categories is > 1, we can use this procedure:
      K <- K / gdc
      init <- rep_len(rep(1:length(K), K), N)
      samples <- unlist(by(init, cats$categories, sample_))
    } else {
      samples <- generate_groups(K)
    }
  } else {
    validate_input(K, "K", objmode = "numeric", len = 1, 
                   greater_than = 1, smaller_than = N)
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



# Generate groups in case of heaviliy uneven group size requirements

# Outer function repeatedly calls generate_groups_() until constraints
# on group sizes are met exactly

generate_groups <- function(tab) {
  loop <- TRUE
  while (loop) {
    samples <- generate_groups_(tab)
    # loop until the number of categories fits the requirements
    loop <- any(sort(table(samples)) != sort(tab))
  }
  samples
}

# Generate group affiliations for all N cases, may not exactly fit 
# the requirements
generate_groups_ <- function(tab) {
  N <- sum(tab) 
  k <- length(tab)
  samples <- unlist(
    replicate(
      N / sum(floor(tab / min(tab))), 
      generate_groups_one_set(tab)
    )
  )
  samples[1:sum(tab)]
}

# Generate one set of group affiliations
generate_groups_one_set <- function(tab) {
  k <- length(tab)
  proportions <- tab / min(tab)
  has_to_be_in <- floor(proportions)
  samples <- rep(1:k, has_to_be_in)
  add <- which(as.logical(rbinom(k, 1, proportions - has_to_be_in)))
  sample(c(samples, add))
}

# Functions to find greatest common denominator among a set of values
gcd <- function(x, y) ifelse(y, Recall(y, x %% y), x)
gcd_ <- function(x) gcd(x[1], x[2])

gcd_all <- function(x) {
  apply(combn(x, 2), 2, gcd_)
}

gcd_set <- function(x) {
  gdc_found <- FALSE 
  while (!gdc_found) {
    x <- gcd_all(x)
    if (all(x == x[1])) {
      gdc_found <- TRUE
    }
  }
  x[1]
}
