
#' Matching
#'
#' Conduct K-partite and unrestricted minimum distance matching to find pairs 
#' or larger groups of similar elements.
#'
#' @param features A numeric vector, matrix or data.frame of data
#'     points.  Rows correspond to elements and columns correspond to
#'     features. A vector represents a single numeric feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. An N x N matrix representing the
#'     pairwise dissimilarities between N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param p The size of the groups; the default is 2, in which case the
#'     function returns pairs.
#' @param groups An optional vector inducing grouping restrictions. If the
#'     vector is passed, the argument \code{p} is ignored and matches are
#'     sought between elements of different groups.
#'
#' @return An integer vector encoding the matches.
#'
#'
#' @details
#' 
#' If the argument \code{features} is passed, matching is done based on the Euclidean
#' distance between data points. To find matches, this function uses
#' the same algorithm as implemented in \code{\link{balanced_clustering}}
#' but it differs with regard to the interface: This function specifies the size of the groups,
#' \code{\link{balanced_clustering}} specifies 
#' the number of the clusters. Moreover, this function makes it possible to specify grouping 
#' restrictions using the \code{groups} argument, thus enabling K-partite minimum matching.
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' # Find triplets
#' N <- 300
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' triplets <- matching(lds, p = 3)
#' plot_clusters(lds, clustering = triplets, within_connection = TRUE)
#'
#'
#' # Match between different plant species
#' species <- iris$Species != "setosa"
#' matched <- matching(iris[species, 1], groups = iris[species, 5])
#' plot_clusters(
#'   data.frame(Species = as.numeric(iris[species, 5]), Sepal.Length = iris[species, 1]),
#'   clustering = matched,
#'   within_connection = TRUE
#' )
#' 
#' @export
#'

matching <- function(features, distances, p = NULL, groups = NULL) {
  data <- process_input(features, distances)
  nn_centroid_clustering(data, p, groups)
}
