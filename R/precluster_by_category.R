
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

# Combining preclustering restrictions and categorical restrictions
# is only possible if
# a) the number of elements per category can be divided by K
# b) the number of elements per anticluster/category combination can be divided by K
preclustering_possible <- function(categories, K) {
  tab <- table(categories)
  !any(tab %% (K^2) != 0)
}
