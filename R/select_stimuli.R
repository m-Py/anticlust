 
# TODO
# - add categorical constrains (using an argument `balance`)
#     + should be able to include preclustering retrictions
# - for subset anticlustering: do categorical preclustering at the
#   beginnung; then, select only preclusters having the correct size;
#   the select the best preclusters
# - add control parameters fro the function select_stimuli (number 
#   of elements per group (as matrix?); p can be control parameter;
#   thresholds on divide parameter; objective for anticlustering 
#   [anticluster editing / k-means]; something to guide importance
#   of maximizing dissimilarity wrt independent variable and minimize
#   similarity wrt covariates)
# - examples in code (need new data set!)
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
#' @param balance Character vector, the names of the variables in 
#'     \code{data} that constitute categorical variables whose frequency
#'     should be balanced across sets
#' @param design Specifies the number of groups per \code{split_by} 
#'     feature. Is a vector of length \code{ncol(split_by)} (or of length
#'     1 if only one - or no - \code{split_by} feature is passed).
#' @param n The number of elements per set.
#' @param randomness An integer in `1:3`. Used for subset anticlustering.
#'    1 = the same stimuli are returned each call; 2 = there is some 
#'    randomness in which stimuli are returned. 3 = there is more 
#'    randomness in which stimuli are returned. 
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
  n = NULL,
  randomness = 1
) {
  message_method(data, split_by, n, design)
  
  if (argument_exists(n)) {
    groups <- subset_anticlustering(data, split_by, equalize, balance, design, n, randomness)
  } else if (argument_exists(split_by) && !argument_exists(n)) {
    groups <- min_max_anticlustering(data, split_by, equalize, balance, design)
  } else if (!argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_anticlustering(data, equalize, balance, design)
  }
  # return the original data with column describing the stimulus set.
  # remove non-selected stimuli
  data$SET <- groups
  data[!is.na(groups), ]
}

# Internal function for subset selection based on preclustering - then 
# anticlustering. if `split_by` exists, min-max anticlustering is conducted
subset_anticlustering <- function(data, split_by, equalize, balance, design, n, randomness) {
  
  N <- nrow(data)
  # track which items are selected via ID:
  data$id_anticlustering <- 1:N

  # get most similar items, possibly under categorical restrictions
  preclusters <- categorical_restrictions(data, equalize, balance, design)
  
  # only select clusters of size `design` (aka `K`)
  filled <- which(table(preclusters) == design)
  is_selected <- preclusters %in% filled
  preselection <- data[is_selected, ]
  preclusters <- preclusters[is_selected]

  # Get best clusters - similar wrt `equalize`; dissimilar wrt `split_by`
  distances <- by(preselection[, equalize], preclusters, dist)
  objectives <- sapply(distances, sum)
    
  # Some additional work needs to be done if min-max anticlustering is required
  if (argument_exists(split_by)) {
    distances2 <- by(preselection[, split_by], preclusters, dist)
    objectives2 <- sapply(distances2, sum)
    objectives <- objectives - objectives2
  }

  needed_n <- design * n
  needed_clusters <- needed_n / table(preclusters)[1]
  # Best preclusters
  if (randomness == 1) {
    cluster_ids <- 1:needed_clusters
  } else if (randomness == 2) {
    # some randomness: better chance to be selected for better objectives
    tmp <- objectives + min(objectives) + 1
    cluster_ids <- sample(length(objectives), size = needed_clusters, prob = sort(tmp / sum(tmp)))
  } else if (randomness == 3) {
    # much randomness: select ANY clusters
    cluster_ids <- sample(length(objectives), size = needed_clusters)
  } else {
    stop("Argument `randomness` must be one of 1, 2, or 3")
  }
  
  most_similar_clusters <- as.numeric(names(sort(objectives))[cluster_ids])
  
  # encode items that are selected (all that are in the best preclusters)
  is_in_output <- preclusters %in% most_similar_clusters
  preselection <- preselection[is_in_output, ]
  preclusters  <- preclusters[is_in_output]

  iv <- NULL
  if (argument_exists(split_by)) {
    iv <- scale(preselection[, split_by])
  }
  
  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    features = scale(preselection[, equalize]),
    K = design,
    iv = iv,
    categories = preclusters
  )
  
  # prepare output: Needs to incorporate NA for non-selected items
  output <- rep(NA, N)
  output[preselection$id_anticlustering] <- anticlusters
  output
}

# A wrapper to obtain preclusters 
# (a) for imbalanced data
# (b) if categorical restrictions are passed
categorical_restrictions <- function(data, equalize, balance, design) {
  if (argument_exists(balance)) {
    categories <- merge_into_one_variable(data[, balance])
    preclusters <- generate_exchange_partners(
      p = design - 1,
      categories = categories,
      similar = TRUE,
      features = scale(data[, equalize])
    )
  } else {
    preclusters <- imbalanced_preclustering(
      scale(data[, equalize]), 
      K = design
    )
  }
  preclusters
}

# Internal function for min-max anticlustering
min_max_anticlustering <- function(data, split_by, equalize, balance, design) {
  K <- prod(design)
  N <- nrow(data)
  preclusters <- categorical_restrictions(data, equalize, balance, design)
  anticlustering(
    features = scale(data[, equalize]),
    K = K,
    iv = scale(data[, split_by]),
    categories = preclusters
  )
}

# Internal function for anticlustering
wrap_anticlustering <- function(data, equalize, balance, design) {
  preclusters <- categorical_restrictions(data, equalize, balance, design)
  
  # initial assignment based on preclustering
  K <- anticlustering(
    data[, equalize],
    K = design,
    categories = preclusters
  )
  
  categories <- NULL
  if (argument_exists(balance)) {
    categories <- merge_into_one_variable(data[, balance])
  }
  
  anticlustering(
    features = scale(data[, equalize]),
    K = K,
    categories = categories
  )
}

message_method <- function(data, split_by, n, design) {
  if (argument_exists(split_by)) {
    min_max <- "Min-Max-"
  } else {
    min_max <- ""
  }
  if (argument_exists(n)) {
    sub <- "Subset-"
  } else {
    sub <- ""
  }
  N <- nrow(data)
  
  if (argument_exists(n)) {
    n <- n
  } else {
    n <- N / design
  }
  if (floor(n) != ceiling(n)) {
    n <- paste(floor(n), "or", ceiling(n))
  }

  message("Starting stimulus selection using `", 
          sub, min_max, "Anticlustering`.")
  message("Selecting ", design, " groups, each having ", 
          n, " elements from a pool of ", N, " stimuli.")
}

