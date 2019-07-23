
#' An objective function measuring similarity of sets
#'
#' Defines set similarity as the discrepancy in means, medians, and
#' standard deviations.
#'
#' @importFrom stats sd var median
#'
#' @param clusters A clustering vector
#' @param features The features
#'
#' @return A value quantifying set similarity. Higher values indicate
#'     that the discrepancy in means, medians and standard deviations
#'     are low.
#'
#' @details
#'
#' Use this function as the \code{objective} in the anticlustering
#' function to make means, standard deviations and medians as similar
#' as possible.
#'
#'
#' @examples
#' data(schaper2019)
#' features <- schaper2019[1:48, 3:6]
#' ac <- anticlustering(
#'   features,
#'   K = 3,
#'   categories = schaper2019$room[1:48],
#'   objective = mean_sd_obj
#' )
#' by(features, ac, function(x) round(colMeans(x), 2))
#' by(features, ac, function(x) round(apply(x, 2, sd), 2))
#' by(features, ac, function(x) round(apply(x, 2, median), 2))
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

## Convert output from `by` to matrix and compute difference based
## on means / sds etc for all features
featurewise_diff <- function(x, K) {
  mat <- matrix(unlist(x), ncol = K)
  mat <- apply(mat, 1, var)
  sum(mat)
}
