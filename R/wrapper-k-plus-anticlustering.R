#' K-plus anticlustering
#' 
#' Perform anticlustering using the k-plus objective to maximize between-group 
#' similarity. This function implements the k-plus anticlustering method described 
#' in Papenberg (2023; <doi:10.1111/bmsp.12315>).
#'
#' @param x A feature matrix where rows correspond to elements and columns
#'     correspond to variables (a single numeric variable can be
#'     passed as a vector).
#' @param K How many anticlusters should be created. Alternatively:
#'     (a) A vector describing the size of each group, or (b) a vector
#'     of length \code{nrow(x)} describing how elements are assigned
#'     to anticlusters before the optimization starts.
#' @param variance Boolean: Should the k-plus objective include a term to 
#'   maximize between-group similarity with regard to the variance? 
#'   (Default = TRUE)
#' @param skew Boolean: Should the k-plus objective include a term to 
#'   maximize between-group similarity with regard to skewness? 
#'   (Default = FALSE)
#' @param kurtosis Boolean: Should the k-plus objective include a term to 
#'   maximize between-group similarity with regard to kurtosis? 
#'   (Default = FALSE)
#' @param covariances Boolean: Should the k-plus objective include a term to 
#'   maximize between-group similarity with regard to covariance structure? 
#'   (Default = FALSE)
#' @param T Optional argument: An integer specifying how many
#'   distribution moments should be equalized between groups. 
#' @param standardize Boolean. If \code{TRUE}, the data is standardized through 
#'     a call to \code{\link{scale}} before the optimization starts. 
#'     Defaults to TRUE. See details.
#' @param ... Arguments passed down to \code{\link{anticlustering}}. All of the 
#'   arguments are supported except for \code{objective}. 
#' 
#' @details 
#'     
#'     This function implements the unweighted sum approach for k-plus 
#'     anticlustering. Details are given in Papenberg (2023). 
#'     
#'     The optional argument \code{T} denotes the number of distribution 
#'     moments that are considered in the anticlustering process. For example,
#'     \code{T = 4} will lead to similar means, variances, skew and kurtosis. 
#'     For the first four moments, it is also possible to use the boolean
#'     convenience arguments \code{variance}, \code{skew} and \code{kurtosis}; the
#'     mean (the first moment) is always included and cannot be "turned off".
#'     If the argument \code{T} is used, it overrides the arguments
#'     \code{variance}, \code{skew} and \code{kurtosis} (corresponding to
#'     the second, third and fourth moment), ignoring their values.
#'     
#'     The \code{standardization} is applied to all original features and the 
#'     additional k-plus features that are appended to the data set in order 
#'     to optimize the k-plus criterion. When using standardization, 
#'     all criteria such as means, variances and skewness receive a comparable
#'     weight during the optimization. It is usually recommended not
#'     to change the default setting \code{standardization = TRUE}.
#'    
#'     This function can use any arguments that are also possible in 
#'     \code{\link{anticlustering}}
#'     (except for `objective` because the objective optimized here
#'     is the k-plus objective; to use a different objective, 
#'     call \code{\link{anticlustering}} directly). Any arguments that are
#'     not explicitly changed here (i.e., \code{standardize = TRUE}) receive the 
#'     default given in \code{\link{anticlustering}} 
#'     (e.g., \code{method = "exchange"}.)
#'     
#' @examples 
#' 
#' # Generate some data
#' N <- 180
#' M <- 4
#' features <- matrix(rnorm(N * M), ncol = M)
#' # standard k-plus anticlustering: optimize similarity with regard to mean and variance:
#' cl <- kplus_anticlustering(features, K = 3, method = "local-maximum")
#' mean_sd_tab(features, cl)
#' # Visualize an anticlustering solution:
#' plot(features, col = palette()[2:4][cl], pch = c(16:18)[cl])
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
#' # The following two calls are equivalent: 
#' init_clusters <- sample(rep_len(1:3, nrow(features)))
#' # 1.
#' x1 <- kplus_anticlustering(
#'   features, 
#'   K = init_clusters, 
#'   variance = TRUE,
#'   skew = TRUE
#' )
#' # 2. 
#' x2 <- kplus_anticlustering(
#'   features, 
#'   K = init_clusters, 
#'   T = 3
#' )
#' # Verify: 
#' all(x1 == x2)
#' 
#' @export
#' 
#' @references
#' 
#' Papenberg, M. (2023). K-plus Anticlustering: An Improved k-means Criterion 
#' for Maximizing Between-Group Similarity. British Journal of Mathematical 
#' and Statistical Psychology. Advance online publication. 
#' https://doi.org/10.1111/bmsp.12315
#' 


kplus_anticlustering <- function(
    x, K, 
    variance = TRUE,
    skew = FALSE,
    kurtosis = FALSE,
    covariances = FALSE,
    T = NULL,
    standardize = TRUE,
    ...) {
  
  validate_input_kplus(x, K, variance, skew, kurtosis, covariances, T, ...)
  
  x <- as.matrix(x)
  M <- ncol(x) 
  feature_bools <- c(variance, skew, kurtosis)
  
  # Unfortunately, I have to extract all ellipsis arguments by hand because
  # preclustering / categories must be used here before calling anticlustering()...
  # (I know this is really ugly...)
  preclustering <- get_argument_from_ellipsis("preclustering", ...)
  preclustering <- ifelse(is.null(preclustering), FALSE, preclustering) # default = FALSE
  categories <- get_argument_from_ellipsis("categories", ...) # default = NULL
  categories <- get_categorical_constraints(x, K, preclustering, categories)
  method <- get_argument_from_ellipsis("method", ...) 
  method <- ifelse(is.null(method), "exchange", method) # default = "exchange"
  repetitions <- get_argument_from_ellipsis("repetitions", ...) # default = NULL

  if (is.null(T)) {
    moments <- (2:4)[feature_bools]
  } else {
    moments <- 2:T
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
    objective = "variance", 
    categories = categories,
    method = method,
    repetitions = repetitions
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

#' Get argument from ellipsis
#' @param argument Name of the argument
#' @param ... The list of additional arguments
#' @examples
#' 
#' get_argument_from_ellipsis(argument = "moep", objective = 1:4)
#' get_argument_from_ellipsis(argument = "objective", objective = 1:4)
#' 
#' @noRd
#' 

get_argument_from_ellipsis <- function(argument, ...) {
  args <- as.list(substitute(list(...)))[-1L]
  if (!argument %in% names(args)) {
    return(NULL)
  }
  eval(args[[argument]])
}


validate_input_kplus <- function(x, K, variance, skew, kurtosis, covariances, T, standardize, ...) {
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
  if (!is.null(T)) {
    validate_input(T, "T", greater_than = 1, must_be_integer = TRUE, len = 1,
                   not_na = TRUE, not_function = TRUE)
  }
  elipsis_arg_names <- names(lapply(substitute(list(...))[-1], deparse))
  if ("objective" %in% elipsis_arg_names) {
    stop("You cannot set the `objective` argument - this is k-plus anticlustering!",
    " \n(-> k-plus is the objective, an extension to k-means anticlustering, see the documenation).")
  }
}
