#' K-plus anticlustering
#' 
#' Perform anticlustering using the k-plus objective to maximize betwee-group 
#' similarity.
#'
#' @param x A feature matrix where rows correspond to elements and columns
#'     correspond to variables (a single numeric variable can be
#'     passed as a vector).
#' @param K How many anticlusters should be created. Alternatively:
#'     (a) A vector describing the size of each group, or (b) a vector
#'     of length \code{nrow(x)} describing how elements are assigned
#'     to anticlusters before the optimization starts.
#' @param variance Boolean: Should the k-plus objective include a term to 
#'   maximizie between-group similarity with regard to the variance? 
#'   (Default = TRUE)
#' @param skew Boolean: Should the k-plus objective include a term to 
#'   maximizie between-group similarity with regard to skewness? 
#'   (Default = FALSE)
#' @param skew Boolean: Should the k-plus objective include a term to 
#'   maximizie between-group similarity with regard to kurtosis? 
#'   (Default = FALSE)
#' @param covariances Boolean: Should the k-plus objective include a term to 
#'   maximizie between-group similarity with regard to covariance structure? 
#'   (Default = FALSE)
#' @param moments An integer vector specifying which distribution moments  
#'   should be equalized between groups.
#' @param standardize Boolean. If \code{TRUE}, the data is standardized through 
#'     a call to \code{\link{scale}} before the optimization starts. 
#'     Defaults to TRUE. See details.
#' @param ... Arguments passed down to \code{link{anticlustering}}. All of the 
#'   arguments are supported except for \code{objective}. Any arguments that are
#'   not explicitly changed here (i.e., \code{standardize}) receive the same 
#'   defaults as in \code{link{anticlustering}}.
#' 
#' @details 
#'     The standardization is applied 
#'     to all original features and the additional features that are appended
#'     in order to optimize the k-plus criterion (this means that all criteria
#'     such as means, variances, skewness etc. receive the same weight during
#'     the optimization.)


k_plus_anticlustering <- function(
    x, K, 
    variance = TRUE,
    skew = FALSE,
    kurtosis = FALSE,
    covariances = FALSE,
    moments = NULL,
    standardize = TRUE
    ...) {
  
}

# function to compute features for variance
# moment = 2 -> variance; 3 = skew; 4 = kurtosis
moment_features <- function(data, moment = 2) {
  apply(data, 2, function(x) (x - mean(x))^moment)
}

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


