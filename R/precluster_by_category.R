
# Compute preclusters within categories
precluster_per_category <- function(features, categories, K) {
  if (preclustering_possible(categories, K)) {
    # save order an restore
    data <- cbind(categories, 1:nrow(features), features)
    data <- sort_by_col(data, 1)
    clusters_by_catgory <- by(data[, -c(1, 2)], data[, 1], balanced_clustering, K = K)
    data$clusters <- unlist(clusters_by_catgory)
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
