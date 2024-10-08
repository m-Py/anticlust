
#' (Anti)cluster editing "diversity" objective
#'
#' Compute the diversity for a given clustering.
#'
#' @param x The data input. Can be one of two structures: (1) A data matrix
#'     where rows correspond to elements and columns correspond to
#'     features (a single numeric feature can be passed as a vector). (2)
#'     An N x N matrix dissimilarity matrix; can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param clusters A vector representing (anti)clusters (e.g.,
#'     returned by \code{\link{anticlustering}}).
#'
#' @return The cluster editing objective
#'
#' @details
#'
#' The objective function used in (anti)cluster editing is the
#' diversity, i.e., the sum of the pairwise distances between elements
#' within the same groups. When the input \code{x} is a feature
#' matrix, the Euclidean distance is computed as the basic distance
#' unit of this objective.
#'
#'
#' @examples
#'
#' data(iris)
#' distances <- dist(iris[1:60, -5])
#' ## Clustering
#' clusters <- balanced_clustering(distances, K = 3)
#' # This is low:
#' diversity_objective(distances, clusters)
#' ## Anticlustering
#' anticlusters <- anticlustering(distances, K = 3)
#' # This is higher:
#' diversity_objective(distances, anticlusters)
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @references
#'
#' Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining
#' diversity and dispersion criteria for anticlustering: A bicriterion
#' approach. British Journal of Mathematical and Statistical
#' Psychology, 73, 275-396. https://doi.org/10.1111/bmsp.12186
#' 
#' Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition 
#' data sets into equivalent parts. Psychological Methods, 26(2), 
#' 161–174. https://doi.org/10.1037/met0000301.
#'
#'

diversity_objective <- function(x, clusters) {
  x <- as.matrix(x)
  validate_input(x, "x", objmode = "numeric")
  validate_input(clusters, "clusters", len = nrow(x), not_na = TRUE)
  diversity_objective_(clusters, x)
}

# other order of the arguments, needed for some internal handling
diversity_objective_ <- function(clusters, x) {
  sum(diversity_objective_by_group(clusters, x))
}

# for local-maximum and multiple repetition method with "average diversity" objective:
weighted_diversity_objective_ <- function(x, clusters, frequencies) {
  sum(diversity_objective_by_group(clusters, x) / frequencies)
}

# Compute distance objective by cluster
# param data: distance matrix or feature matrix
# param cl: cluster assignment
diversity_objective_by_group <- function(cl, data) {
  if (is_distance_matrix(data)) {
    objectives <- sapply(
      sort(unique(cl)), 
      function(x) sum(as.dist(data[cl == x, cl == x]))
    )
  } else {
    objectives <- sapply(
      sort(unique(cl)), 
      function(x) sum(dist(data[cl == x, , drop = FALSE]))
    )
  }
  objectives
}
