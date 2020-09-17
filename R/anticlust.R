#' anticlust: Subset Partitioning via Anticlustering
#'
#' The method of anticlustering partitions a pool of elements into
#' groups (i.e., anticlusters) in such a way that the between-group
#' similarity is maximized and -- at the same time -- the within-group
#' heterogeneity is maximized. This reverses the logic of cluster
#' analysis that strives for high within-group homogeneity and low
#' similarity of the different groups. Computationally, anticlustering
#' is accomplished by maximizing instead of minimizing a clustering
#' objective function, such as the intra-cluster variance (used in
#' k-means clustering) or the sum of pairwise distances within
#' clusters.  The function anticlustering() implements exact and
#' heuristic anticlustering algorithms as described in Papenberg and
#' Klau (2020; <doi:10.1037/met0000301>). The exact approach requires
#' that the GNU linear programming kit
#' (<https://www.gnu.org/software/glpk/glpk.html>) is available and
#' the R package 'Rglpk' (<https://cran.R-project.org/package=Rglpk>)
#' is installed. Some other functions are available to solve classical
#' clustering problems. The function balanced_clustering() applies a
#' cluster analysis under size constraints, i.e., creates equal-sized
#' clusters. The function matching() can be used for (unrestricted,
#' bipartite, or K-partite) matching. The function wce() can be used
#' optimally solve the (weighted) cluster editing problem, also known
#' as correlation clustering, clique partitioning problem or
#' transitivity clustering.
#' 
#' @section Primary functions:
#' \code{\link{anticlustering}} 
#' \code{\link{balanced_clustering}} 
#' \code{\link{matching}} 
#' \code{\link{categorical_sampling}} 
#' \code{\link{wce}}
#'
#' @docType package
#' @name anticlust
#' @useDynLib anticlust, .registration = TRUE
#' 
#' 
NULL
