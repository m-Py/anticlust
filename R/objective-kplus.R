
#' Objective value for the k-plus criterion
#'
#'
#' @param x A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param clusters A vector representing (anti)clusters (e.g.,
#'     returned by \code{\link{anticlustering}} or
#'     \code{\link{balanced_clustering}})
#'
#' @return The value of the k-plus criterion.
#'
#' @details
#' 
#' The k-plus criterion is an extension of the k-means criterion
#' (i.e., the "variance", see \code{\link{variance_objective}}). The
#' standard k-means objective is high if the means of the input
#' variables are similar between clusters, but there is no guarantee
#' that the standard deviations will also be similar (in fact,
#' maximizing the k-means objective tends to decrease similarity in
#' standard deviations in comparison to a completely random
#' assignment). However, to achieve overall between-group similarity,
#' it is desirable that the spread of the data is also similar between
#' groups---and not just the means. This is accomplished by maximizing
#' the k-plus criterion that also incorporates the standard deviations
#' of the input variables.
#'
#' Equalizing means and standard deviations simultaneously is
#' accomplished by internally appending new variables to the data
#' input \code{x}, one new variable for each column in \code{x}. These
#' new variables contain the squared difference of each data point to
#' the mean of the respective column, and are then included---in
#' addition to the original data---in standard k-means
#' anticlustering. This way, the average squared deviation of the data
#' points to the means becomes similar between groups, which is the
#' variance. Hence, the k-plus criterion simultaneously represents
#' similarity in means and variance (and thus, the standard
#' deviation), and can be used to simultaneously equalize the mean and
#' the spread of the data.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @note
#'
#' K-plus anticlustering has newly been implemented in the package
#' anticlust (available since version 0.5.2). The author is currently
#' working on a paper detailing the objective's background, but has
#' already made the methodology available as its results have been
#' very convincing thus far (e.g., check out the examples below). When
#' using k-plus anticlustering in your research, it would be courteous
#' to cite Papenberg and Klau (2020) as the primary \code{anticlust}
#' reference, even though the criterion has not been described in that
#' paper. In doubt, contact the author to inquire whether a new
#' reference is available, or check out the package website
#' (https://github.com/m-Py/anticlust).
#'
#' @references
#' 
#' Papenberg, M., & Klau, G. W. (2020). Using anticlustering to partition 
#' data sets into equivalent parts. Psychological Methods. Advance Online 
#' Publication. https://doi.org/10.1037/met0000301.
#' 
#'
#' @examples
#'
#' data(schaper2019)
#' features <- schaper2019[, 3:6]
#' 
#' # Optimize k-plus criterion
#' kplus_groups <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "kplus"
#' )
#' # Optimize normal k-means criterion
#' kmeans_groups <- anticlustering(
#'   features,
#'   K = 3,
#'   objective = "variance"
#' )
#' 
#' # Compute k-plus criterion (k-plus is much better here)
#' kplus_objective(features, kplus_groups)
#' kplus_objective(features, kmeans_groups)
#' 
#' # Compare to k-means criterion (k-plus not much worse here)
#' variance_objective(features, kplus_groups)
#' variance_objective(features, kmeans_groups)
#'
#' # Compare means and standard deviations after k-means and k-plus
#' # anticlustering (the standard deviations are usually much closer
#' # after k-plus anticlustering, but there is only little to no
#' # difference with regard to the means)
#' mean_sd_tab(features, kplus_groups)
#' mean_sd_tab(features, kmeans_groups)
#'
#' 
#' 

kplus_objective <- function(x, clusters) {
  x <- as.matrix(x)
  variance_objective(cbind(x, squared_from_mean(x)), clusters)
}

squared_from_mean <- function(data) {
  apply(data, 2, function(x) (x - mean(x))^2)
}
