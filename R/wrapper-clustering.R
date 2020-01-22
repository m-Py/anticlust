
#' Create balanced clusters of equal size
#'
#' @param x The data input. Can be one of two structures: (1) A data matrix
#'     where rows correspond to elements and columns correspond to
#'     features (a single numeric feature can be passed as a vector). (2)
#'     An N x N matrix dissimilarity matrix; can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param K How many clusters should be created.
#' @param method One of "heuristic" or "ilp". See details.
#'
#' @return An integer vector representing the cluster affiliation of 
#'     each data point
#'
#' @details
#'
#' This function partitions a set of elements into \code{K} equal-sized
#' clusters. The function offers two methods, a heuristic method and an
#' exact method. The heuristic (\code{method = "heuristic"}) computes
#' the centroid of all available elements and identifies the element
#' farthest to it. If the input is a dissimilarity matrix, the most
#' central element acts as the centroid. The farthest element is
#' clustered with its \code{(N/K) - 1} nearest neighbours. From the remaining
#' elements, the element farthest to the centroid is selected and again
#' clustered with its \code{(N/K) - 1} neighbours; the procedure is repeated
#' until all elements are part of a cluster. 
#'
#' An exact method (\code{method = "ilp"}) can be used to solve cluster
#' editing optimally. The cluster editing objective minimizes the sum
#' of pairwise distances within clusters. If the argument
#' \code{features} is passed, the Euclidean distance is computed by
#' default as the basic unit of the cluster editing objective. If
#' another distance measure is preferred, users may pass a self-computed
#' dissimiliarity matrix via the argument \code{distances}. The optimal
#' cluster editing objective can be found via integer linear
#' programming. To obtain an optimal solution for balanced cluster
#' editing, a linear programming solver must be installed and usable
#' from R. The `anticlust` package supports the open source GNU linear
#' programming kit (called from the package \code{Rglpk}) and the
#' commercial solvers gurobi (called from the package \code{gurobi}) and
#' IBM CPLEX (called from the package \code{Rcplex}). A license is
#' needed to use one of the commercial solvers.
#'
#'
#' @source
#'
#' The heuristic method was originally developed and contributed by m.eik michalke.
#' It was later rewritten by Martin Papenberg, who also implemented the exact integer linear
#' programming method.
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' m.eik michalke \email{meik.michalke@@hhu.de}
#'
#' @examples
#'
#' # Cluster a data set and visualize results
#' N <- 1000
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' cl <- balanced_clustering(lds, K = 10)
#' plot_clusters(lds, clusters = cl)
#' 
#' # Repeat using a distance matrix as input
#' cl2 <- balanced_clustering(dist(lds), K = 10)
#' plot_clusters(lds, clusters = cl2)
#'
#' @references
#'
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59–96.
#'

balanced_clustering <- function(x, K, method = "heuristic") {

  input_validation_anticlustering(x, K, "distance", method, TRUE, NULL)
  
  data <- process_input(x)
  
  if (method == "ilp") {
    return(balanced_cluster_editing(data, K, solver_available()))
  }
  nn_centroid_clustering(data, NROW(data) / K)
}
