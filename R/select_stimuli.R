 
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
# - control parameter: add some randomness to increase diversity of 
#   output (this may be important for subset anticlustering)
 
#' Select stimuli for experiments
#' 
#' Stimulus selection via Anticlustering or Min-Max-Anticlustering
#' 
#' @param data A N x M data frame of features describing stimuli
#' @param split_by Character vector, the names of the variables in 
#'     \code{data} that should be different between sets
#' @param equalize Character vector, the names of the variables in 
#'     \code{data} that should be similar across sets
#' @param balance One or more categorical variables. The frequency of
#'     the categories is balanced between sets.
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

select_stimuli <- function(
  data, 
  split_by = NULL, 
  equalize, 
  balance = NULL,
  design, 
  n = NULL
) {
  if (argument_exists(split_by) && argument_exists(n)) {
    groups <- subset_minmax_anticlustering(data, split_by, equalize, design, n)
  } else if (argument_exists(split_by) && !argument_exists(n)) {
    groups <- min_max_anticlustering(data, split_by, equalize, design)
  } else if (!argument_exists(split_by) && argument_exists(n)) {
    groups <- subset_anticlustering(data, equalize, balance, design, n)
  } else if (!argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_anticlustering(data, equalize, design)
  }
  # return the original data with column describing the stimulus set.
  # remove non-selected stimuli
  data$SET <- groups
  data[!is.na(groups), ]
}

# Internal function for subset selection based on preclustering - then anticlustering
subset_minmax_anticlustering <- function(data, split_by, equalize, design, n) {
  N <- nrow(data)
  message("Starting stimulus selection using `Subset Min-Max-Anticlustering`.")
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
  
  # Get best clusters - similar wrt `equalize`; dissimilar wrt `split_by`
  distances <- by(preselection[, equalize], preclusters, dist)
  distances2 <- by(preselection[, split_by], preclusters, dist)
  objectives <- sapply(distances, sum)
  objectives2 <- sapply(distances2, sum)
  objectives <- objectives - objectives2
  
  needed_n <- design * n
  needed_clusters <- needed_n / table(preclusters)[1]
  # Best preclusters
  most_similar_clusters <- as.numeric(names(sort(objectives))[1:needed_clusters])
  
  # encode items that are selected (all that are in the best preclusters)
  is_in_output <- preclusters %in% most_similar_clusters
  preselection <- preselection[is_in_output, ]
  preclusters  <- preclusters[is_in_output]

  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    features = scale(preselection[, equalize]),
    K = prod(design),
    iv = scale(preselection[, split_by]),
    categories = preclusters
  )
  
  # prepare output: Needs to incorporate NA for non-selected items
  output <- rep(NA, N)
  output[preselection$id_anticlustering] <- anticlusters
  output
}

# Internal function for subset selection based on preclustering - then anticlustering
subset_anticlustering <- function(data, equalize, balance, design, n) {
  N <- nrow(data)
  message("Starting stimulus selection using `Subset Anticlustering`.")
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

  # if a categorical variable is passed, use this for the initial 
  # assignment and exchange method; can't use both categorical and
  # preclustering restrictions (only in very specific cases - I should
  # implement the test `preclustering_possible()` and 
  # `precluster_per_category()` for this case
  if (argument_exists(balance)) {
    preclusters  <- merge_into_one_variable(preselection[, balance])
  }
  
  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    scale(preselection[, equalize]),
    K = design,
    categories = preclusters
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
  message("Starting stimulus selection using `Min-Max Anticlustering`.")
  message("Selecting ", K, " groups, each having approximately ", 
          N / K, " elements from a pool of ", N, " stimuli.")
  
  preclusters <- imbalanced_preclustering(scale(data[, equalize]), nrow(data) / K)
  
  anticlustering(
    features = scale(data[, equalize]),
    K = K,
    iv = scale(data[, split_by]),
    categories = preclusters
  )
}

# Internal function for anticlustering
wrap_anticlustering <- function(data, equalize, design) {
  message("Starting stimulus selection using `Anticlustering`.")
  message("Selecting ", design, " groups, each having approximately ", 
          round(nrow(data) / design), " elements from a pool of ", 
          nrow(data), " stimuli.")
  
  # initial assignment based on preclustering
  K <- anticlustering(
    data[, equalize],
    K = design,
    categories = imbalanced_preclustering(scale(data[, equalize]), nrow(data) / design)
  )
  
  anticlustering(
    features = scale(data[, equalize]),
    K = K
  )
}
