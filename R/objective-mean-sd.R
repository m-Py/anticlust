
#' An objective function measuring similarity of sets
#'
#' Compute the discrepancy in means and standard deviations between clusters.
#'
#' @importFrom stats sd var
#'
#' @param features A matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features.
#' @param clusters A clustering vector
#'
#' @return A value quantifying similarity in means and standard deviations. 
#'     Higher values indicate that means and standard deviations are more
#'     similar.
#'
#' @details
#'
#' This function can be passed as the argument \code{objective} to the 
#' function \code{\link{anticlustering}} to minimize differences in means 
#' and standard deviations between anticlusters.
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
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#'

mean_sd_obj <- function(features, clusters) {
   K <- length(unique(clusters))
   mean_min  <- featurewise_diff(by(features, clusters, colMeans, na.rm = TRUE), K)
   SD_min <- featurewise_diff(by(features, clusters, function(x) apply(x, 2, var, na.rm = TRUE)), K)
   (mean_min + SD_min) * (-1)
}

## Convert output from `by` to matrix and compute difference based
## on means / sds etc for all features
featurewise_diff <- function(x, K) {
  mat <- matrix(unlist(x), ncol = K)
  mat <- apply(mat, 1, var)
  sum(mat)
}
