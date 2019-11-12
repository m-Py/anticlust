 
# TODO
# - use preclustering in anticlustering applications 
#   (generate_exchange_partners with `similar = TRUE`)
# - add control parameters fro the function select_stimuli (number 
#   of elements per group (as matrix?); p can be control parameter;
#   thresholds on divide parameter; objective for anticlustering 
#   [anticluster editing / k-means]; something to guide importance
#   of maximizing dissimilarity wrt independent variable and minimize
#   similarity wrt covariates)
# - examples in code (need new data set!)
# - argument p by default should behave differently for divide and select
#   and anticlustering; anticlustering: no restriction on exchange 
#   partners. divide and select: 15 exchange partners
# - generate_exchange_partners should call the clustering function defined
#   below that can deal with data that cannot be split in N/K parts
# - control parameter for weights of `equalize` variables
# - maybe: argument pca: should a pca be conducted on the features?
# - use categories argument (as argument `balance`)
# - 
 
#' Select stimuli for experiments
#' 
#' Stimulus selection via `Divide and Select` or `anticlustering`.
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
#' @param p The number of exchange partners per stimulus used by the
#'     exchange method that optimizes similarity with regard to the 
#'     \code{equalize} variables. Higher values increase
#'     the precision of the results but also increase run time.
#'
#' @return The grouping of each item
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

select_stimuli <- function(
  data, 
  split_by = NULL, 
  equalize, 
  balance,
  design, 
  n = NULL
) {
  if (argument_exists(split_by) && argument_exists(n)) {
    groups <- divide_and_select(data, split_by, equalize, design, n)
  } else if (argument_exists(split_by) && !argument_exists(n)) {
    groups <- min_max_anticlustering(data, split_by, equalize, design)
  } else if (!argument_exists(split_by) && argument_exists(n)) {
    groups <- subset_anticlustering(data, equalize, design, n)
  } else if (!argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_anticlustering(data, equalize, design)
  }
  groups
}

# Internal function for subset selection based on preclustering - then anticlustering
subset_anticlustering <- function(data, equalize, design, n) {
  N <- nrow(data)
  message("Starting stimulus selection using `preclustering and anticlustering`.")
  message("Selecting ", design, " groups, each having ", 
          n, " elements from a pool of ", N, " stimuli.")

  # keep track which items are selected; use ID for this
  data$id_anticlustering <- 1:N
  # Select a number of elements that can evenly be split into preclusters
  preselection <- data[sample(1:N, size = N - (N %% design)), ]

  # Compute preclusters
  preclusters <- balanced_clustering(
    scale(preselection[, equalize]), 
    K = nrow(preselection) / design
  )
  
  # Get most similar clusters
  distances <- by(preselection[, equalize], preclusters, dist)
  objectives <- sapply(distances, sum)
  needed_n <- design * n
  needed_clusters <- needed_n / table(preclusters)[1]
  # Best preclusters
  most_similar_clusters <- as.numeric(names(sort(objectives))[1:needed_clusters])
  
  # encode items that are selected (all that are in the best preclusters)
  is_in_output <- preclusters %in% most_similar_clusters
  preselection <- preselection[is_in_output, ]
  preclusters  <- preclusters[is_in_output]
  
  # First assignment: balance preclusters across anticlusters!
  K <- anticlustering(
    preselection[, equalize], 
    K = design, 
    method = "sampling", 
    nrep = 1,
    categories = preclusters
  )
  
  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    scale(preselection[, equalize]),
    K = K
  )
  
  # prepare output: Needs to incorporate NA for non-selected items
  output <- rep(NA, N)
  output[preselection$id_anticlustering] <- anticlusters
  output
}

# Internal function for min-max anticlustering
min_max_anticlustering <- function(data, split_by, equalize, design) {
  K <- prod(design)
  N <- nrow(data)
  message("Starting stimulus selection using `min-max anticlustering`.")
  message("Selecting ", K, " groups, each having approximately ", 
          N / K, " elements from a pool of ", N, " stimuli.")
  anticlustering(
    features = scale(data[, equalize]),
    K = K,
    iv = scale(data[, split_by])
  )
}

# Internal function for anticlustering
wrap_anticlustering <- function(data, equalize, design) {
  message("Starting stimulus selection using `anticlustering`.")
  message("Selecting ", design, " groups, each having approximately ", 
          round(nrow(data) / design), " elements from a pool of ", 
          nrow(data), " stimuli.")
  anticlustering(
    features = scale(data[, equalize]),
    K = design
  )
}

# Internal function for divide and select approach
divide_and_select <- function(data, split_by, equalize, design, n) {
  message("Starting stimulus selection using the `Divide and Select` approach.")
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
