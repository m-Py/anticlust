 
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
  validate_input_selection(data, split_by, equalize, balance, design, n) 
  message_method(data, split_by, n, design)
  if (argument_exists(n)) {
    groups <- subset_anticlustering(data, split_by, equalize, balance, design, n)
  } else {
    preclusters <- preclustering(data, equalize, balance, prod(design))
    groups <- wrap_anticlustering(data, equalize, split_by, design, preclusters)
  }
  # return the original data with column describing the stimulus set.
  # remove non-selected stimuli
  data$SET <- groups
  data[!is.na(groups), ]
}

wrap_anticlustering <- function(data, equalize, split_by, design, preclusters) {
  # First, generate objective function (either anticlustering or 
  # min-max-anticlustering); default objective: anticluster editing
  obj_fun <- "distance"
  # min-max anticlustering if `split_by` exists
  if (argument_exists(split_by)) {
    obj_fun <- make_obj_function(data, equalize, split_by, design) 
  }
  
  # optimize set assignment using exchange method  
  anticlusters <- anticlustering(
    features = scale(data[, c(equalize, split_by)]),
    K = prod(design),
    categories = preclusters,
    objective = obj_fun
  )
}


# Generate an objective functions for min-max anticlustering
make_obj_function <- function(data, equalize, split_by, design) {
  
  # get all levels of clusters and the corresponding levels for split variables
  levels <- expand.grid(lapply(design, function(x) 1:x))
  
  # combine to a single objective function
  function(cl, data) {
    similarity_covariates <- obj_value_distance(
      cl, 
      data[, equalize, drop = FALSE]
    )
    dissims <- c()
    for (i in 1:length(design)) {
      clusters <- levels[, i][cl]
      # compute sum of distances within cluster
      distances <- by(data[, split_by[i]], clusters, dist)
      dissims[i] <- sum(sapply(distances, sum))
    }
    similarity_covariates - sum(dissims)
  }
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

