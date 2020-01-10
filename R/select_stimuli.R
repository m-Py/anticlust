 
# TODO
# - generate_exchange_partners should call the clustering function defined
#   below that can deal with data that cannot be split in N/K parts
# - control parameter for weights of `equalize` / `split_by` variables


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
subset_anticlustering <- function(data, split_by, equalize, balance, design, n) {
  
  # total N of input data:
  N <- nrow(data)
  # number of groups that are needed:
  K <- prod(design)
  
  # track which items are selected via ID:
  data$id_anticlustering <- 1:N
  
  # test if I can safely remove all entries with NA values
  preselection <- safely_exclude_na(data, split_by, equalize, balance, K, n)

  # get most similar items, possibly under categorical restrictions
  preclusters <- categorical_restrictions(preselection, equalize, balance, K)
  
  # only select clusters of size `K`
  filled <- which(table(preclusters) == K)
  is_selected <- preclusters %in% filled
  preselection <- preselection[is_selected, , drop = FALSE]
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

  needed_n <- K * n
  needed_clusters <- needed_n / table(preclusters)[1]
  # Best preclusters
  cluster_ids <- 1:needed_clusters
  most_similar_clusters <- as.numeric(names(sort(objectives))[cluster_ids])
  
  # encode items that are selected (all that are in the best preclusters)
  is_in_output <- preclusters %in% most_similar_clusters
  preselection <- preselection[is_in_output, ]
  preclusters  <- preclusters[is_in_output]

  iv <- NULL
  obj_fun <- "distance"
  if (argument_exists(split_by)) {
    iv <- scale(preselection[, split_by])
    obj_fun <- make_obj_function(data, equalize, split_by, design) 
  }
  
  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    features = scale(preselection[, c(equalize, split_by)]),
    K = K,
    categories = preclusters,
    objective = obj_fun
  )
  
  # prepare output: Needs to incorporate NA for non-selected items
  output <- rep(NA, N)
  output[preselection$id_anticlustering] <- anticlusters
  output
}

safely_exclude_na <- function(data, split_by, equalize, balance, K, n) {
  has_no_na <- complete.cases(data[, equalize, drop = FALSE])
  if (argument_exists(split_by)) {
    has_no_na <- has_no_na & complete.cases(data[, split_by, drop = FALSE])
  }
  if (argument_exists(balance)) {
    has_no_na <- has_no_na & complete.cases(data[, balance, drop = FALSE])
  }
  # can remove data if enough entries remain after NA exclusion
  if ((K * n) <= sum(has_no_na)) {
    data <- data[has_no_na, , drop = FALSE]
    if (sum(!has_no_na) > 0) {
      message("\n", sum(!has_no_na), " records could be excluded due to missing values.\n",
              "Selecting stimuli from the remaining ", nrow(data), " records.")
    }
  }
  data
}

# A wrapper to obtain preclusters 
# (a) for imbalanced data
# (b) if categorical restrictions are passed
categorical_restrictions <- function(data, equalize, balance, K) {
  if (argument_exists(balance)) {
    categories <- merge_into_one_variable(data[, balance])
    preclusters <- generate_exchange_partners(
      p = K - 1,
      categories = categories,
      similar = TRUE,
      features = scale(data[, equalize])
    )
  } else {
    preclusters <- imbalanced_preclustering(
      scale(data[, equalize]), 
      K = K
    )
  }
  preclusters
}

# Internal function for min-max anticlustering
min_max_anticlustering <- function(data, split_by, equalize, balance, design) {
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

