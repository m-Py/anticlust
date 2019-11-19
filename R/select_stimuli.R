 
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
  design, 
  n = NULL
) {
  if (argument_exists(split_by) && argument_exists(n)) {
    groups <- subset_anticlustering(data, split_by, equalize, design, n)
  } else if (argument_exists(split_by) && !argument_exists(n)) {
    groups <- min_max_anticlustering(data, split_by, equalize, design)
  } else if (!argument_exists(split_by) && argument_exists(n)) {
    groups <- subset_anticlustering(data, split_by, equalize, design, n)
  } else if (!argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_anticlustering(data, equalize, design)
  }
  # return the original data with column describing the stimulus set.
  # remove non-selected stimuli
  data$SET <- groups
  data[!is.na(groups), ]
}

# Internal function for subset selection based on preclustering - then 
# anticlustering. if `split_by` exists, min-max anticlustering is conducted
subset_anticlustering <- function(data, split_by, equalize, design, n) {
  N <- nrow(data)

  # keep track which items are selected; use ID for this
  data$id_anticlustering <- 1:N
  # Select a number of elements that can evenly be split into preclusters
  preselection <- remove_outliers(data, equalize, prod(design))

  # Compute preclusters
  preclusters <- balanced_clustering(
    scale(preselection[, equalize]), 
    K = nrow(preselection) / design
  )
  
  # Get best clusters - similar wrt `equalize`; dissimilar wrt `split_by`
  distances <- by(preselection[, equalize], preclusters, dist)
  objectives <- sapply(distances, sum)
    
  # Some additional work needs to be done if min-max anticlustering is required
  min_max <- ""
  if (argument_exists(split_by)) {
    min_max <- "Min-Max-"
    distances2 <- by(preselection[, split_by], preclusters, dist)
    objectives2 <- sapply(distances2, sum)
    objectives <- objectives - objectives2
  }

  needed_n <- design * n
  needed_clusters <- needed_n / table(preclusters)[1]
  # Best preclusters - (TODO: insert stochastic component)
  most_similar_clusters <- as.numeric(names(sort(objectives))[1:needed_clusters])
  
  # encode items that are selected (all that are in the best preclusters)
  is_in_output <- preclusters %in% most_similar_clusters
  preselection <- preselection[is_in_output, ]
  preclusters  <- preclusters[is_in_output]

  iv <- NULL
  if (argument_exists(split_by)) {
    iv <- scale(preselection[, split_by])
  }

  message("Starting stimulus selection using `Subset-", min_max, "Anticlustering`.")
  message("Selecting ", design, " groups, each having ", 
          n, " elements from a pool of ", N, " stimuli.")
  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    features = scale(preselection[, equalize]),
    K = prod(design),
    iv = iv,
    categories = preclusters
  )
  
  # prepare output: Needs to incorporate NA for non-selected items
  output <- rep(NA, N)
  output[preselection$id_anticlustering] <- anticlusters
  output
}

# Function to remove elements such that the remaining elements can fit 
# into clusters of size K. Removes the data points that are furthest 
# apart from any other data points
remove_outliers <- function(data, equalize, K) {
  N <- nrow(data)
  if (N %% K == 0) {
    return(data)
  }
  distances <- as.matrix(dist(data[, equalize]))
  diag(distances) <- Inf
  minima <- apply(distances, 1, min)
  minima <- cbind(minima, 1:N)
  minima <- sort_by_col(minima, 1)
  # discard elements whose nearest neighbor is farthest away
  minima <- minima[1:(N - (N %% K)), , drop = FALSE]
  # obtain original order
  minima <- sort_by_col(minima, 2)
  data[minima[, 2], ]
}

# Internal function for min-max anticlustering
min_max_anticlustering <- function(data, split_by, equalize, design) {
  K <- prod(design)
  N <- nrow(data)
  message("Starting stimulus selection using `Min-Max Anticlustering`.")
  message("Selecting ", K, " groups, each having approximately ", 
          N / K, " elements from a pool of ", N, " stimuli.")
  
  preclusters <- imbalanced_preclustering(scale(data[, equalize]), K)
  
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
    categories = imbalanced_preclustering(scale(data[, equalize]), design)
  )
  
  anticlustering(
    features = scale(data[, equalize]),
    K = K
  )
}
