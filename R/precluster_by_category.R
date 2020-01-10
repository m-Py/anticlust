
# A wrapper to obtain preclusters 
# (a) for imbalanced data
# (b) if categorical restrictions are passed
categorical_restrictions <- function(data, equalize, balance, K) {
  if (argument_exists(balance)) {
    categories <- merge_into_one_variable(data[, balance])
    preclusters <- to_numeric(
      precluster_per_category(
        data[, equalize, drop = FALSE], 
        categories, 
        K
      )
    )
  } else {
    preclusters <- imbalanced_preclustering(
      scale(data[, equalize]), 
      K = K
    )
  }
  preclusters
}

# Compute preclusters within categories
# argument K is the number of anticlusters, not the number of clusters!
precluster_per_category <- function(features, categories, K) {
  N <- nrow(features)
  unique_categories <- sort(unique(categories))
  n_categories <- length(unique_categories)
  # save original order to restore before returning
  data <- cbind(categories, 1:N, features)
  data <- sort_by_col(data, 1)
  data <- data.frame(data)
  # precluster within each category:
  cl <- list()
  for (i in 1:n_categories) {
    tmp_data <- data[data[, 1] == unique_categories[i], -c(1, 2), drop = FALSE]
    rownames(tmp_data) <- NULL
    if (nrow(tmp_data) <= K) {
      cl[[i]] <- 1 # all members have the same category and cluster
    } else {
      cl[[i]] <- imbalanced_preclustering(tmp_data, K = K)
    }
  }
  # ensure that clusters in different categories have different cluster numbers
  to_be_added <- sapply(cl, length) / K
  to_be_added <- cumsum(to_be_added) - to_be_added[1]
  cl <- lapply(seq_along(cl), function(i) cl[[i]] + to_be_added[i])
  # return data sorted in original order
  data$clusters <- unlist(cl)
  cl <- sort_by_col(data, 2)$clusters
  merge_into_one_variable(cbind(cl, categories))
}

# Wrap the balanced clustering function to allow for some imbalancing
# @param features
# @param K the Number of anticlusters (! -> the number of clusters returned is ~ nrow(features) / K)

imbalanced_preclustering <- function(features, K) {
  features <- data.frame(features)
  N <- nrow(features)
  features$preclusters_id <- 1:N
  preclusters <- rep(NA, N)
  # only select as many data as can be clustered into balanced clusters
  subsetted <- remove_outliers(features, K)
  # test if only one cluster can be filled:
  if (nrow(subsetted) == K)  {
    preclusters[subsetted$preclusters_id] <- 1
  } else {
    cl <- nn_centroid_clustering(
      subsetted[, -ncol(subsetted)],
      K = K
    )
    preclusters[subsetted$preclusters_id] <- cl
  }
  # the remainders get into the last cluster
  if (sum(is.na(preclusters)) > 0) {
    preclusters[is.na(preclusters)] <- max(preclusters, na.rm = TRUE) + 1
  }
  preclusters
}


# Function to remove elements such that the remaining elements can fit 
# into clusters of size K. Removes the data points that are furthest 
# apart from the data centroid
remove_outliers <- function(data, K) {
  N <- nrow(data)
  if (N %% K == 0) {
    return(data)
  }
  # define outliers as being furthest away from data centroid
  distances <- distances_from_centroid(data)
  # get indices of elements that are no outliers
  idx <- order(distances)[1:(N - (N %% K))]
  # return these elements
  data[sort(idx), , drop = FALSE]
}
