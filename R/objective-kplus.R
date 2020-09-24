
#' Objective value for the k-plus criterion
#'
#'
#' @param x A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param clusters A vector representing (anti)clusters (e.g., returned
#'     by \code{\link{anticlustering}} or
#'     \code{\link{balanced_clustering}})
#'
#' @return The value of the k-plus criterion.
#'
#' @details
#' 
#' The k-plus criterion is an extension of the k-means objective (i.e., "variance")
#' (see \code{\link{variance_objective}}). The standard k-means objective is high
#' the means of the input variables are similar between groups, but there is no
#' guarantee that the standard deviations will also be similar (in fact, optimizing
#' the k-means objective alone tends to decrease similarity in standard deviations). 
#' The k-plus criterion
#' also tries to make the standard deviations similar between groups. This is 
#' accomplished by internally appending new variables to the data input \code{x}, 
#' one new variable for each column in \code{x}. These new variables contain the squared 
#' difference of each data point to the center of the respective column. If 
#' these variables are included in standard k-means anticlustering, 
#' their means tend to become similar between groups. This means that the
#' mean squared deviation of the data points to 
#' the mean becomes similar between groups, i.e., the variance. Hence,
#' the k-plus criterion simultaneously represents similarity in means and 
#' variance / standard deviation and can be used to simultaneously optimize 
#' both objectives.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' data(schaper2019)
#' features <- schaper2019[, 3:6]
#' 
#' # Optimize k-plus criterion
#' groups1 <- anticlustering(
#'   features,
#'   K = 4,
#'   objective = "kplus"
#' )
#' # Optimize normal k-means criterion
#' groups2 <- anticlustering(
#'   features,
#'   K = 4,
#'   objective = "variance"
#' )
#' 
#' # Compute k-plus criterion (k-plus is much better here)
#' kplus_objective(features, groups1)
#' kplus_objective(features, groups2)
#' 
#' # Compare to k-means criterion (k-plus not much worse here)
#' variance_objective(features, groups1)
#' variance_objective(features, groups2)
#' 

kplus_objective <- function(x, clusters) {
  x <- as.matrix(x)
  variance_objective(cbind(x, squared_from_mean(x)), clusters)
}

squared_from_mean <- function(data) {
  apply(data, 2, function(x) (x - mean(x))^2)
}
