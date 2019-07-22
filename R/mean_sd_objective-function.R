
#' An objective function measuring similarity of sets
#'
#' Uses variance in means, medians, and standard deviations per feature
#' across sets as the objective.
#'
#' @importFrom stats sd var median
#'
#' @param clusters A clustering vector
#' @param features The features
#'
#' @return A value quantifying set similarity as the total absolute deviations
#'     between sets on all features. Higher values indicate that the deviations
#'     in mean, median and standard deviations are low.
#'
#' @export
#'

mean_sd_obj <- function(clusters, features) {
   K <- length(unique(clusters))
   mean_min  <- featurewise_diff(by(features, clusters, colMeans), K)
   median_min  <- featurewise_diff(by(features, clusters, function(x) apply(x, 2, median)), K)
   SD_min <- featurewise_diff(by(features, clusters, function(x) apply(x, 2, sd)), K)
   (mean_min + SD_min + median_min) * (-1)
}

## Convert output from `by` to matrix and compute difference based on means / sds etc
## for all features
featurewise_diff <- function(x, K) {
  mat <- matrix(unlist(x), ncol = K)
  mat <- apply(mat, 1, var)
  sum(mat)
}
