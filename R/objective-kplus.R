
#' Compute k-plus variables
#'
#'
#' @param x A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param T The number of distribution moments for which variables are generated.
#' @param standardize Logical, should all columns of the output be standardized
#'  (defaults to TRUE).
#'
#' @return A matrix containing all columns of \code{x} and all additional
#' columns of k-plus variables. If \code{x} has M columns, the output matrix
#' has M * T columns.
#'
#' @details
#' 
#' The k-plus criterion is an extension of the k-means criterion
#' (i.e., the "variance", see \code{\link{variance_objective}}).
#' In \code{\link{kplus_anticlustering}}, equalizing means and variances
#' simultaneously (and possibly additional distribution moments) is
#' accomplished by internally appending new variables to the data
#' input \code{x}. When using only the  variance as additional criterion, the
#' new variables represent the squared difference of each data point to
#' the mean of the respective column. All columns are then included---in
#' addition to the original data---in standard k-means
#' anticlustering. The logic is readily extended towards higher order moments,
#' see Papenberg (2024). This function gives users the possibility to generate
#' k-plus variables themselves, which offers some additional flexibility when
#' conducting k-plus anticlustering.
#'
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @references
#' 
#' Papenberg, M. (2024). K-plus Anticlustering: An Improved k-means Criterion for 
#' Maximizing Between-Group Similarity. British Journal of Mathematical and 
#' Statistical Psychology, 77(1), 80--102. https://doi.org/10.1111/bmsp.12315
#' 
#' @export
#'
#' @examples
#' 
#' # Use Schaper data set for example
#' data(schaper2019)
#' features <- schaper2019[, 3:6]
#' K <- 3
#' N <- nrow(features)
#'
#' # Some equivalent ways of doing k-plus anticlustering:
#' 
#' init_groups <- sample(rep_len(1:3, N))
#' table(init_groups)
#' 
#' kplus_groups1 <- anticlustering(
#'   features,
#'   K = init_groups,
#'   objective = "kplus",
#'   standardize = TRUE,
#'   method = "local-maximum"
#' )
#'
#' kplus_groups2 <- anticlustering(
#'   kplus_moment_variables(features, T = 2), # standardization included by default
#'   K = init_groups,
#'   objective = "variance", # (!)
#'   method = "local-maximum"
#' )
#' 
#' # this function uses standardization by default unlike anticlustering():
#' kplus_groups3 <- kplus_anticlustering(
#'   features, 
#'   K = init_groups,
#'   method = "local-maximum"
#' )
#'
#' all(kplus_groups1 == kplus_groups2)
#' all(kplus_groups1 == kplus_groups3)
#' all(kplus_groups2 == kplus_groups3)
#' 

kplus_moment_variables <- function(x, T, standardize = TRUE) {
  validate_data_matrix(x)
  validate_input(T, "T", greater_than = 1, must_be_integer = TRUE, len = 1,
                 not_na = TRUE, not_function = TRUE)
  validate_input(standardize, "standardize", objmode = "logical", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  x <- as.matrix(x)
  M <- ncol(x)
  for (i in 2:T) {
    x <- cbind(x, moment_features(x, T))
  }
  colnames(x) <- paste0(colnames(x), "^", rep(1:T, each = M))
  if (!standardize) {
    return(x)
  }
  scale(x)
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

squared_from_mean <- function(data) {
  apply(data, 2, function(x) (x - mean(x))^2)
}
