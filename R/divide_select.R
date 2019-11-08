 
#' Divide and Select 
#' 
#' Stimulus selection through the divide and select approach. 
#' Makes some variables dissimilar between sets, other variables
#' similar between sets. 
#'
#' @param split_by The features by which the data is split
#' @param select The features that should be similar across sets
#' @param design Specifies the number of levels per \code{split_by} 
#'     feature. Is a vector of length \code{ncol(split_by)} (or of length
#'     1 if only one \code{split_by} feature is passed.
#' @param n The number of elements per set.
#' @param p The number of exchange partners; higher values increase
#'     the precision of the results but also increase run time.
#'
#' @return The anticluster affiliation
#'
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @export
#' 
#' @examples
#' # TODO
#' 
#' TODO control parameters: 
#'   - number of elements per group (as matrix?)
#'   - p can be control parameter
#'   - thresholds on divide parameter
 

divide_and_select <- function(split_by, equalize, design, n, p) {
  if (argument_exists(split_by)) {
    split_by <- split_data(split_by, design)
  }
  categories <- merge_into_one_variable(split_by)
  init_groups <- initialize_K(groups = categories, n = rep(n, length(unique(categories))))
  exchange_partners <- generate_exchange_partners(categories = categories, p = p)
  anticlustering(
    scale(equalize),
    K = init_groups,
    categories = exchange_partners,
    objective = mean_sd_obj
  )
}

# Split data / Turn numeric values into categories
# 
# param: split_by Inherited from divide_and_select
# param: design Inherited from divide_and_select
# return: a matrix of same dimensions as `split_by` encoding the split per variable
#
split_data <- function(split_by, design) {
  # this function should be able to incorporate custom thresholds passed by the user
  # now it will simply use cutoffs "in the middle"
  split_by <- as.list(data.frame(split_by))
  mapply(categorize_vector, x = split_by, k = design)
}


# Turn a numeric vector into ordinal values (categories)
# param x: the vector
# param k: The number of categories
# return: the categorized vector
categorize_vector <- function(x, k) {
  x <- matrix(c(1:length(x), x), ncol = 2)
  x <- sort_by_col(x, 2)
  x <- cbind(x, sort(rep_len(1:k, nrow(x))))
  # return categorized data in original order
  sort_by_col(x, 1)[, 3]
}

