
#' Matching
#'
#' Conduct K-partite or unrestricted (minimum distance) matching to find pairs 
#' or groups of similar elements. Finding matches is by default based on the Euclidean 
#' distance between data points, but a custom distance measure can also be employed.
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
#' N <- 299
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' triplets <- matching(lds, p = 3)
#' plot_clusters(lds, clustering = triplets, within_connection = TRUE)
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
#' # Use unequal-sized groups: Only selects matches for some elements
#' N <- 33
#' data <- matrix(rnorm(N), ncol = 1)
#' groups <- sample(1:2, size = N, replace = TRUE, prob = c(0.7, 0.3))
#' matched <- matching(data[, 1], groups = groups)
#' plot_clusters(
#'   cbind(groups, data), 
#'   clustering = matched, 
#'   within_connection = TRUE
#' )
#' 
#' 
#' @export
#'

matching <- function(features = NULL, distances = NULL, p = 2, groups = NULL) {
  data <- process_input(features, distances)
  N <- nrow(data)
  groups <- merge_into_one_variable(groups)
  cl <- nn_centroid_clustering(data, p, groups)
  # Before returning: order the group numbers by objective - most similar 
  # matches have lower indices
  sort_by_objective(cl, data, N)
}

sort_by_objective <- function(cl, data, N) {
  selected <- (1:N)[!is.na(cl)]
  cl_sub <- cl[selected]
  cl_sub <- order_cluster_vector(cl_sub)
  if (is_distance_matrix(data)) {
    data <- data[selected, selected]
    objectives <- sapply(
      1:max(cl_sub), 
      function(x) sum(as.dist(data[cl_sub == x, cl_sub == x]))
    )
  } else {
    data <- data[selected, , drop = FALSE]
    objectives <- sapply(
      1:max(cl_sub), 
      function(x) sum(dist(data[cl_sub == x, ]))
    )
  }
  # recode original matching labels according to objective
  N <- length(cl_sub)
  # sorry for this code, but it works and it was really complicated
  one <- data.frame(match = cl_sub, order_matches = 1:length(cl_sub))
  two <- data.frame(match = order(objectives), objective_id = 1:length(objectives))
  merged <- merge(one, two)
  new <- rep(NA, N)
  new[merged$order_matches] <- merged$objective_id
  cl[!is.na(cl)] <- new
  cl
}
