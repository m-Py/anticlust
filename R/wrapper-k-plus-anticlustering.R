#' K-plus anticlustering
#' 
#' Perform anticlustering using the k-plus objective to maximize between-group 
#' similarity. This function implements the k-plus anticlustering method described 
#' Papenberg (2023; <doi:10.31234/osf.io/7jw6v>).
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
#' @param kurtosis Boolean: Should the k-plus objective include a term to 
#'   maximizie between-group similarity with regard to kurtosis? 
#'   (Default = FALSE)
#' @param covariances Boolean: Should the k-plus objective include a term to 
#'   maximizie between-group similarity with regard to covariance structure? 
#'   (Default = FALSE)
#' @param moments Optional argument: An integer vector specifying which 
#'   distribution moments should be equalized between groups.
#' @param standardize Boolean. If \code{TRUE}, the data is standardized through 
#'     a call to \code{\link{scale}} before the optimization starts. 
#'     Defaults to TRUE. See details.
#' @param ... Arguments passed down to \code{\link{anticlustering}}. All of the 
#'   arguments are supported except for \code{objective}. 
#' 
#' @details 
#'     
#'     If the argument \code{moments} is used, it overrides the arguments
#'     \code{variance}, \code{skew} and \code{kurtosis} (corresponding to
#'     the second, third and fourth moment), ignoring their values. Note that 
#'     the first moment, i.e., the mean is always included in the optimization
#'     (corresponding to the standard k-means criterion); it cannot be "turned off"
#'     by using the argument \code{moments}.
#'     
#'     The \code{standardization} is applied to all original features and the 
#'     additional features that are appended in order to optimize 
#'     the k-plus criterion (this means that all criteria such as means, 
#'     variances, skewness etc. receive the same relative weight 
#'     during the optimization.)
#'    
#'     This function can use any arguments that are also possible in 
#'     \code{\link{anticlustering}}
#'     (except for `objective` of course; because the objective optimized here
#'     is the k-plus objective -- to use a different objective, 
#'     call \code{\link{anticlustering}} directly). Any arguments that are
#'     not explicitly changed here (i.e., \code{standardize = TRUE}) receive the 
#'     default given in \code{\link{anticlustering}} 
#'     (e.g., `method = "exchange"`.)
#'     
#' @examples 
#' 
#' # Generate some data
#' N <- 180
#' M <- 4
#' features <- matrix(rnorm(N * M), ncol = M)
#' plot(features)
#' # standard k-plus anticlustering: optimize similarity with regard to mean
#' # and variance:
#' cl <- kplus_anticlustering(features, K = 3, method = "local-maximum")
#' mean_sd_tab(features, cl)
#' 
#' # Also optimize with regard to skewness and kurtosis
#' cl2 <- kplus_anticlustering(
#'   features, 
#'   K = 3, 
#'   method = "local-maximum", 
#'   skew = TRUE, 
#'   kurtosis = TRUE
#' )
#' 
#' # Try to equalize the first 10 moments between groups (the first moment, 
#' # i.e., the mean, is always considered in k-plus anticlustering)
#' kplus_anticlustering(
#'   features, 
#'   K = 3, 
#'   moments = 2:10, 
#' )
#' @export
#' 
#' @references
#' 
#' Papenberg, M. (2023). k-plus Anticlustering: An Improved k-means Criterion 
#' for Maximizing Between-Group Similarity. Retrieved from psyarxiv.com/7jw6v 
#' 


kplus_anticlustering <- function(
    x, K, 
    variance = TRUE,
    skew = FALSE,
    kurtosis = FALSE,
    covariances = FALSE,
    moments = NULL,
    standardize = TRUE,
    ...) {
  
  validate_input_kplus(x, K, variance, skew, kurtosis, covariances, moments, ...)
  
  x <- as.matrix(x)
  M <- ncol(x) 
  feature_bools <- c(variance, skew, kurtosis)
  

  if (is.null(moments)) {
    moments <- (2:4)[feature_bools]
  }
  # determine the number of features after augmentation, for validation below
  M_augmented <- M + M * length(moments) + ifelse(covariances, choose(M, 2), 0)

  # Add additional k-plus features: 
  #   1. Moments
  augmented_data <- x
  for (i in seq_along(moments)) {
    augmented_data <- cbind(augmented_data, moment_features(x, moments[i]))
  }
  #   2. Covariance    
  if (covariances == TRUE) {
    augmented_data <- cbind(augmented_data, covariance_features(x))
  }
  
  # Validation of number of features
  stopifnot(ncol(augmented_data) == M_augmented)
  
  # Call anticlustering
  anticlustering(
    augmented_data, 
    K = K, 
    standardize = standardize, 
    objective = "variance", ...
  )

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


validate_input_kplus <- function(x, K, variance, skew, kurtosis, covariances, moments, standardize, ...) {
  validate_data_matrix(x)
  # K is validated in anticlustering()
  validate_input(variance, "variance", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  validate_input(skew, "skew", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  validate_input(kurtosis, "kurtosis", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  validate_input(covariances, "covariances", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  validate_input(moments, "moments", greater_than = 1, must_be_integer = TRUE,
                 not_na = TRUE, not_function = TRUE)
  elipsis_arg_names <- names(lapply(substitute(list(...))[-1], deparse))
  if ("objective" %in% elipsis_arg_names) {
    stop("You cannot set the `objective` argument - this is k-plus anticlustering!",
    " \n(-> k-plus is the objective, an extension to k-means anticlustering, see the documenation).")
  }
}
