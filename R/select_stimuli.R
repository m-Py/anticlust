 
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
#'     1 if one or no \code{split_by} feature is passed).
#' @param n The number of elements per set. If not specified, the complete 
#'     date set is partitioned into subsets (the number of subsets is specified
#'     using the \code{design} argument).
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
  if (argument_exists(split_by)) {
    if (length(design) != length(split_by)) {
      stop("Length of argument `design` must match length of argument `split_by`.")
    }
  }
  message_method(data, split_by, n, design)
  if (argument_exists(n)) {
    groups <- subset_anticlustering(data, split_by, equalize, balance, design, n)
  } else if (argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_min_max_anticlustering(data, split_by, equalize, balance, design)
  } else if (!argument_exists(split_by) && !argument_exists(n)) {
    groups <- wrap_anticlustering(data, equalize, balance, design)
  }
  # return the original data with column describing the stimulus set.
  # remove non-selected stimuli
  data$SET <- groups
  data[!is.na(groups), ]
}

# Internal function for min-max anticlustering
wrap_min_max_anticlustering <- function(data, split_by, equalize, balance, design) {
  K <- prod(design)
  N <- nrow(data)
  preclusters <- categorical_restrictions(data, equalize, balance, K)
  anticlustering(
    features = scale(data[, equalize]),
    K = K,
    iv = scale(data[, split_by]),
    categories = preclusters
  )
}

# Internal function for anticlustering
wrap_anticlustering <- function(data, equalize, balance, K) {
  preclusters <- categorical_restrictions(data, equalize, balance, K)
  
  # initial assignment based on preclustering
  K <- anticlustering(
    data[, equalize],
    K = K,
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
  K <- prod(design)
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
    n <- N / K
  }
  if (floor(n) != ceiling(n)) {
    n <- paste(floor(n), "or", ceiling(n))
  }

  message("Starting stimulus selection using `", 
          sub, min_max, "Anticlustering`.")
  message("Selecting ", K, " groups, each having ", 
          n, " elements from a pool of ", N, " stimuli.")
}

