 
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
# - argument p b default should behave differently for divide and select
#   and anticlustering; anticlustering: no restriction on exchange 
#   partners. divide and select: 15 exchange partners
# - generate_exchange_partners should call the clustering function defined
#   below that can deal with data that cannot be split in N/K parts
 
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

select_stimuli <- function(data, split_by = NULL, equalize, design, n = NULL, p = 15) {
  if (argument_exists(split_by) && argument_exists(n)) {
    groups <- divide_and_select(data, split_by, equalize, design, n, p)
  } else if (argument_exists(split_by) && !argument_exists(n)) {
    groups <- min_max_anticlustering(data, split_by, equalize, design, p)
  } else if (!argument_exists(split_by) && argument_exists(n)) {
    groups <- subset_anticlustering(data, equalize, design, n, p)
  } else if (!argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_anticlustering(data, equalize, design, p)
  }
  groups
}

# Internal function for subset selection based on preclustering - then anticlustering
subset_anticlustering <- function(data, equalize, design, n, p) {
  N <- nrow(data)
  message("Starting stimulus selection using `preclustering and anticlustering`.")
  message("Selecting ", design, " groups, each having ", 
          n, " elements from a pool of ", N, " stimuli.")

  # create vector that keeps track which item is in the selection
  preclusters <- rep(NA, N)
  mode(preclusters) <- "numeric"
  gets_precluster <- rep(TRUE, N)
  # randomly discard some items so that a balanced selection is possible
  preselection <- sample(1:N, size = N - (N %% design))
  not_selected <- (1:N)[!(1:N) %in% preselection]
  gets_precluster[not_selected] <- FALSE
  
  # compute preclusters
  subsetted <- data[preselection, ]
  N <- nrow(subsetted)
  preclusters[gets_precluster] <- balanced_clustering(scale(subsetted[, equalize]), K = N / design)
  # compute sum of within-cluster distances per cluster
  distances <- by(subsetted[, equalize], preclusters[!is.na(preclusters)], dist)
  objective <- sapply(distances, sum)
  # how many elements remain?
  needed_n <- design * n
  # how many clusters are needed:
  needed_clusters <- needed_n / table(preclusters)[1]
  # select the best clusters
  most_similar_clusters <- as.numeric(names(sort(objective))[1:needed_clusters])
  is_not_in_output <- which(!preclusters %in% most_similar_clusters)
  is_in_output <- which(preclusters %in% most_similar_clusters)
  preclusters[is_not_in_output] <- NA
  
  # now do anticlustering
  anticlusters <- rep(NA, nrow(data))
  subsetted <- data[is_in_output, ]
  preclusters <- preclusters[is_in_output]
  
  # First assignment: balance preclusters across anticlusters!
  K <- anticlustering(
    subsetted[, equalize], 
    K = design, 
    method = "sampling", 
    nrep = 1,
    categories = preclusters
  )
  
  # now optimize assignment using exchange method
  anticlusters[is_in_output] <- anticlustering(
    scale(subsetted[, equalize]),
    K = K, 
    categories = generate_exchange_partners(nrow(subsetted), p = p)
  )
  anticlusters
}

# Internal function for min-max anticlustering
min_max_anticlustering <- function(data, split_by, equalize, design, p) {
  K <- prod(design)
  N <- nrow(data)
  message("Starting stimulus selection using `min-max anticlustering`.")
  message("Selecting ", K, " groups, each having approximately ", 
          N / K, " elements from a pool of ", N, " stimuli.")
  # Get preclusters based on p argument. 
  k <- p + 1
  preclusters <- imbalanced_preclustering(data, N, k, equalize)
  anticlustering(
    features = scale(data[, equalize]),
    K = K,
    categories = preclusters,
    iv = data[, split_by]
  )
}

# Deal with it that it may not be possible to generate perfectly 
# balanced clusters
imbalanced_preclustering <- function(data, N, k, equalize) {
  preclusters <- rep(NA, N)
  # only select as many data as can be clustered into balanced clusters
  subsetted <- data[1:(N - (N %% k)), ]
  preclusters[1:nrow(subsetted)] <- balanced_clustering(
    scale(subsetted[, equalize]), 
    K = nrow(subsetted) / k
  )
  # full some clusters at random
  if (sum(is.na(preclusters)) > 0) {
    preclusters[is.na(preclusters)] <- sample(1:max(preclusters, na.rm = TRUE), size = sum(is.na(preclusters)))
  }
  preclusters
}

# Internal function for anticlustering
wrap_anticlustering <- function(data, equalize, design, p) {
  message("Starting stimulus selection using `anticlustering`.")
  message("Selecting ", design, " groups, each having approximately ", 
          round(nrow(data) / design), " elements from a pool of ", 
          nrow(data), " stimuli.")
  anticlustering(
    features = scale(data[, equalize]),
    K = design,
    categories = generate_exchange_partners(nrow(data), p = p)
  )
}

# Internal function for divide and select approach
divide_and_select <- function(data, split_by, equalize, design, n, p) {
  message("Starting stimulus selection using the `Divide and Select` approach.")
  message("Selecting ", prod(design), " groups, each having ", n, 
          " elements from a pool of ", nrow(data), " stimuli.")
  equalize <- data[, equalize]
  split_by <- data[, split_by]
  split_by <- split_data(split_by, design)
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
  x <- matrix(c(1:length(x), x), ncol = 2)
  x <- sort_by_col(x, 2)
  x <- cbind(x, sort(rep_len(1:k, nrow(x))))
  # return categorized data in original order
  sort_by_col(x, 1)[, 3]
}

