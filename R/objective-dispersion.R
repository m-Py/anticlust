#' Cluster dispersion 
#' 
#' Compute the dispersion objective for a given clustering (i.e., the
#' minimum distance between two elements within the same cluster).
#' 
#' @param x The data input. Can be one of two structures: (1) A
#'     feature matrix where rows correspond to elements and columns
#'     correspond to variables (a single numeric variable can be
#'     passed as a vector). (2) An N x N matrix dissimilarity matrix;
#'     can be an object of class \code{dist} (e.g., returned by
#'     \code{\link{dist}} or \code{\link{as.dist}}) or a \code{matrix}
#'     where the entries of the upper and lower triangular matrix
#'     represent pairwise dissimilarities.
#' @param clusters A vector representing (anti)clusters (e.g.,
#'     returned by \code{\link{anticlustering}}).
#'     
#' @details 
#' 
#' The dispersion is the minimum distance between two elements within
#' the same cluster. When the input \code{x} is a feature matrix, the
#' Euclidean distance is used as the distance unit. Maximizing the
#' dispersion maximizes the minimum heterogeneity within clusters and
#' is an anticlustering task.
#' 
#' @examples 
#' 
#' N <- 50 # number of elements
#' M <- 2  # number of variables per element
#' K <- 5  # number of clusters
#' random_data <- matrix(rnorm(N * M), ncol = M)
#' random_clusters <- sample(rep_len(1:K, N))
#' dispersion_objective(random_data, random_clusters)
#' 
#' # Maximize the dispersion 
#' optimized_clusters <- anticlustering(
#'   random_data,
#'   K = random_clusters, 
#'   objective = dispersion_objective
#' )
#' dispersion_objective(random_data, optimized_clusters)
#' 
#' @export
#'
#' @references
#'
#' Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining
#' diversity and dispersion criteria for anticlustering: A bicriterion
#' approach. British Journal of Mathematical and Statistical
#' Psychology, 73, 275-396. https://doi.org/10.1111/bmsp.12186

#' 

dispersion_objective <- function(x, clusters) {
  x <- convert_to_distances(x)
  validate_input(x, "x", objmode = "numeric")
  validate_input(clusters, "clusters", len = nrow(x), not_na = TRUE)
  dispersion_objective_(clusters, x)
}

dispersion_objective_ <- function(clusters, distances) {
  objectives <- sapply(
    sort(unique(clusters)), 
    # min() gives warning for clusters of length 0 (but its behaviour is as intended)
    function(x) suppressWarnings(min(as.dist(distances[clusters == x, clusters == x])))
  )
  min(objectives)
}
