
#' An objective function measuring similarity of sets
#'
#' Uses variance in means and standard deviations per feature across sets
#' as the objective.
#'
#' @importFrom stats sd
#'
#' @param clusters A clustering vector
#' @param features The features
#'
#' @return A value quantifying set similarity as the total absolute deviations
#'     between sets on all features.
#'
#' @export
#'

mean_sd_obj <- function(clusters, features) {
   features <- features[!is.na(clusters), , drop = FALSE]
   clusters <- clusters[!is.na(clusters)]
   K <- length(unique(clusters))
   O_min  <- featurewise_diff(by(features, clusters, colMeans), K)
   SD_min <- featurewise_diff(by(features, clusters, function(x) apply(x, 2, sd)), K)
   (O_min + SD_min) * (-1)
}

## Convert output from `by` to matrix and compute difference based on means / sds etc
## for all features
featurewise_diff <- function(x, K) {
  mat <- matrix(unlist(x), ncol = K)
  mat <- apply(mat, 1, diff)
  sum(abs(mat))
}
