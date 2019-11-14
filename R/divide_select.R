
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
  anticlustering(
    scale(equalize),
    K = init_groups,
    categories = exchange_partners,
    objective = mean_sd_obj
  )
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
