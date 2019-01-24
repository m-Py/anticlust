
# Test that the number of features is a multiplier of unique
# anticlusters
legal_number_of_clusters <- function(features, clusters) {
  ## 1. correct number of clusters assignments?
  if (length(clusters) != nrow(features))
    stop("The number of cluster assignments and the number of features have to match.")
  n_anticlusters <- length(unique(clusters))
  ## 2. More than 1 disctinct cluster?
  if (n_anticlusters <= 1)
    stop("There have to be at least two different anticlusters")
  ## 3. Do all clusters occur equally often?
  if (!all(table(clusters) == table(clusters)[1]))
    stop("Each clusters must occur equally often")
  ## 4. Probably redundant to 3:
  if (nrow(features) %% n_anticlusters != 0)
    stop("The number of elements is not a multiplier of the number of anticlusters")
  invisible(NULL)
}
