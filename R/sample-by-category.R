
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

  # If groups have unequal size:
  if (sum(K) == N) {
    gdc <- gcd_set(K)
    K <- K / gdc
    init <- rep_len(rep(1:length(K), K), N)
  } else { # same size for each group
    validate_input(K, "K", objmode = "numeric", len = 1, 
                   greater_than = 1, must_be_integer = TRUE, not_na = TRUE)
    init <- rep_len(1:K, N)
  }
  
  cats <- data.frame(
    categories = categories,
    order = 1:N
  )
  cats <- sort_by_col(cats, "categories")
  cats$samples <- unlist(by(init, cats$categories, sample_))
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
