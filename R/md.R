
#' Compute mahalanobis distance between two vectors
#'
#' @param x First vector
#' @param y First vector
#' @param cov_mat covariance matrix between the features in x and y
#' @param lambda Regularization, defaults to 0 (no regularization)
#'
#' @importFrom MASS ginv
#'
#' @return The mahalanobis distance (scalar value)
#'
#' @export
#'
#' @examples
#' mahalanobis_distance(iris[1, -5], iris[2, -5], cov(iris[, -5]))
#'

mahalanobis_distance <- function(x, y, cov_mat, lambda = 0) {
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == nrow(cov_mat))
  stopifnot(length(x) == ncol(cov_mat))
  M <- length(x)
  MD <- double(1)
  INVERTED_COV_MATRIX <- MASS::ginv(cov_mat+lambda*diag(nrow(cov_mat)))
  results <- .C(
      "mahalanobis_out",
      as.double(x),
      as.double(y),
      as.integer(M),
      as.double(INVERTED_COV_MATRIX),
      MD = as.double(MD),
      PACKAGE = "anticlust"
  )
  results$MD
}
