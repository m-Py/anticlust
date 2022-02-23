
# Compute k-covariances features for an original number of M features.
# That is choose(M, 2) additional features, i.e., O(M^2) additional features.

covariance_features <- function(data) {
  M <- ncol(data)
  if (M < 2) {
    stop("k-covariances can only be used when at least two features are present.")
  }
  feature_combinations <- t(combn(M, 2))
  # store new features in matrix
  cov_feature_matrix <- matrix(ncol = nrow(feature_combinations), nrow = nrow(data))
  for (i in 1:nrow(feature_combinations)) {
    f1 <- feature_combinations[i, 1]
    f2 <- feature_combinations[i, 2]
    cov_feature_matrix[, i] <- (data[, f1] - mean(data[, f1])) * (data[, f2] - mean(data[, f2]))
  }
  cov_feature_matrix
}


