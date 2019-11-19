
#' Select stimuli for experiments
#' 
#' Stimulus selection via `Divide and Select`
#' 
#' @param data A N x M data frame of features describing stimuli
#' @param split_by Character vector, the names of the variables in 
#'     \code{data} that should be different between sets
#' @param equalize Character vector, the names of the variables in 
#'     \code{data} that should be similar across sets
#' @param design Specifies the number of groups per \code{split_by} 
#'     feature. Is a vector of length \code{ncol(split_by)} (or of length
#'     1 if only one - or no - \code{split_by} feature is passed).
#' @param n The number of elements per set.
#'
#' @return A data frame that has the same columns as the original input 
#'    (the data frame \code{data}), but has an additional called \code{SET}.
#'
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @export
#' 
#' @examples
#' # TODO
#' 
#' @details
#' The argument \code{split_by} will convert a numeric variable into 
#' a categorical variable, and will recognize if a variable is already
#' categorical by testing if \code{length(unique(split_by)) == design}
#' is \code{TRUE}.
#' 

# Internal function for divide and select approach
divide_and_select <- function(data, split_by, equalize, design, n) {
  message("Starting stimulus selection using `Divide and Select`.")
  message("Selecting ", prod(design), " groups, each having ", n, 
          " elements from a pool of ", nrow(data), " stimuli.")
  equalize <- data[, equalize]
  split_by <- data[, split_by]
  split_by <- split_data(split_by, design)
  categories <- merge_into_one_variable(split_by)
  init_groups <- initialize_K(groups = categories, n = rep(n, length(unique(categories))))
  exchange_partners <- generate_exchange_partners(categories = categories, p = 15)
  groups <- anticlustering(
    scale(equalize),
    K = init_groups,
    categories = exchange_partners,
    objective = mean_sd_obj
  )
  data$SET <- groups
  data[!is.na(groups), ]
}

# Split data / Turn numeric values into categories
# 
# param: split_by A matrix with variables to be different between sets
# param: design Inherited from `select_stimuli`
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
  # test if this is already a categorical vector
  if (length(unique(x)) == k) {
    return(x)
  }
  x <- matrix(c(1:length(x), x), ncol = 2)
  x <- sort_by_col(x, 2)
  x <- cbind(x, sort(rep_len(1:k, nrow(x))))
  # return categorized data in original order
  sort_by_col(x, 1)[, 3]
}
