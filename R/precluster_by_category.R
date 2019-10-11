
# Compute preclusters within categories
precluster_per_category <- function(features, categories, K) {
  if (preclustering_possible(categories, K)) {
    N <- nrow(features)
    unique_categories <- unique(categories)
    n_categories <- length(unique_categories)
    # save original order to restore before returning
    data <- cbind(categories, 1:N, features)
    data <- sort_by_col(data, 1)
    # precluster within each category:
    cl <- list()
    for (i in 1:n_categories) {
      tmp_data <- data[data[, 1] == unique_categories[i], -c(1, 2)]
      cl[[i]] <- balanced_clustering(tmp_data, K = nrow(tmp_data) / K)
    }
    # ensure that clusters in different categories have different cluster numbers
    to_be_added <- sapply(cl, length) / K
    to_be_added <- cumsum(to_be_added) - to_be_added[1]
    cl <- lapply(seq_along(cl), function(i) cl[[i]] + to_be_added[i])
    # return data sorted in original order
    data$clusters <- unlist(cl)
    return(sort_by_col(data, 2)$clusters)
  }
  stop("Not today.")
}

# Combining preclustering restrictions and categorical restrictions
# is only possible if
# a) the number of elements per category can be divided by K
# b) the number of elements per anticluster/category combination can be divided by K
preclustering_possible <- function(categories, K) {
  tab <- table(categories)
  !any(tab %% (K^2) != 0)
}
