
#' Weighted cluster editing / Clique partitioning
#'
#' Optimally solves the weighted cluster editing problem, also known as
#' »correlation clustering« or »clique partitioning problem«.
#'
#' @param x A N x N similarity matrix. Larger values indicate stronger
#'     agreement / similarity between a pair of data points
#' @param solver Optional argument; if passed, has to be either "glpk"
#'     or "symphony". See details.
#' @param method Either "ilp" (default) or "local-maximum".
#' @param repetitions Number of repetitions when using \code{method =
#'     "local-maximum"}.
#'
#' @return An integer vector representing the cluster affiliation of
#'     each data point
#'
#' @examples
#' \donttest{
#' features <- swiss
#' distances <- dist(scale(swiss))
#' hist(distances)
#' # Define agreement as being close enough to each other according to 
#' # the Euclidean distance. By coding low agreement as -1 and high agreement 
#' # as +1, we solve *unweighted* cluster editing
#' agreements <- ifelse(as.matrix(distances) < 3, 1, -1)
#' clusters <- wce(agreements)
#' plot(swiss, col = clusters, pch = 19)
#' }
#' 
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @details 
#' 
#' Finds the clustering that maximizes the sum of pairwise
#' similarities within clusters. In the input \code{x}, some
#' similarities should be negative (indicating dissimilarity) because
#' otherwise the maximum sum of similarities is obtained by simply
#' joining all elements within a single big cluster. The cluster
#' editing / clique partitioning method has the advantage that the
#' number of clusters does not need to be specified by the user
#' (unlike, for example, in k-means or p-medoids clustering) and
#' instead is found by the algorithm. However, the solution depends on
#' the classification of pairwise similarity versus dissimilarity,
#' indicated by positive and negative values in the similarity matrix
#' \code{x}. Different ways of obtaining a suitable similarity matrix
#' are available (e.g., Brusco & Köhn, 2009; Wittkop et al., 2010).
#' 
#' By default, an optimal ILP method is used, implementing the model
#' by Grötschel and Wakabayashi (1989). If the argument \code{solver}
#' is not specified, the function will try to find the GLPK or
#' SYMPHONY solver by itself (it prioritizes using SYMPHONY if both
#' are available).
#' 
#' The \code{method = "local-maximum"} implements an adaptation of the
#' local neighbourhood search that Späth (1986) proposed for
#' maximizing the k-means criterion. However, the initial clustering
#' is not chosen randomly, but instead each element starts as a
#' separate "singleton" cluster. After that, each element is
#' sequentially assigned all possible cluster labels (that is,
#' \code{nrow(x)} cluster labels), and for each element, the best
#' exchange is realized with regard to the increasing the
#' within-cluster sum of similarities (see
#' \code{\link{diversity_objective}}). The search terminates when no
#' improvement can occur by exchanging the cluster label of an
#' element. When specifying the argument \code{repetitions}, the local
#' maximum search is repeated, but the order in which the elements are
#' processed is changed randomly each time. Note that this algorithm
#' is not deterministic - even when the \code{repetitions} argument is
#' not specified - because the order in which elements are processed is
#' always chosen randomly.
#' 
#' @note
#' 
#' When using \code{method = "ilp"}, this function either requires the
#' R package \code{Rglpk} and the GNU linear programming kit
#' (<http://www.gnu.org/software/glpk/>) or the R package
#' \code{Rsymphony} and the COIN-OR SYMPHONY solver libraries
#' (<https://github.com/coin-or/SYMPHONY>).
#' 
#' @references
#'
#' 
#' Bansal, N., Blum, A., & Chawla, S. (2004). Correlation clustering. 
#' Machine Learning, 56, 89–113. 
#' 
#' Brusco, M. J., & Köhn, H. F. (2009). Clustering qualitative data
#' based on binary equivalence relations: neighborhood search
#' heuristics for the clique partitioning problem. Psychometrika, 74,
#' 685–703.
#' 
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59-96.
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213-218.
#' 
#' Wittkop, T., Emig, D., Lange, S., Rahmann, S., Albrecht, M.,
#' Morris, J. H., ..., Baumbach, J. (2010). Partitioning biological
#' data with transitivity clustering. Nature Methods, 7, 419–420.
#'

wce <- function(x, solver = NULL, method = "ilp", repetitions = NULL) {

  validate_data_matrix(x)
  if (!is_distance_matrix(x)) {
    stop("The input via argument `weights` is not a similarity matrix, ",
         "the upper and lower triangulars of your matrix differ.")
  }
  validate_input(
    method, "method", input_set = c("ilp", "local-maximum"), 
    not_na = TRUE, not_function = TRUE, len = 1, objmode = "character"
  )
  
  x <- as.matrix(x)
  if (method == "local-maximum") {
    if (argument_exists(repetitions)) {
      validate_input(repetitions, "repetitions", objmode = "numeric", len = 1, 
                     greater_than = 0, must_be_integer = TRUE, not_na = TRUE,
                     not_function = TRUE)
    }
    repetitions <- ifelse(is.null(repetitions), 1, repetitions)
    if (repetitions > 1) {
      solutions <- lapply(1:repetitions, heuristic_wce, x = x)
      # Get best of all solutions
      objs <- lapply(
        solutions,
        diversity_objective_,
        x = x
      )
      return(solutions[[which.max(objs)]])
    } else {
      return(heuristic_wce(NULL, x))
    }
  } 
  # exact cluster editing:
  if (argument_exists(solver)) {
    validate_input(solver, "solver", objmode = "character", len = 1,
                   input_set = c("glpk", "symphony"), not_na = TRUE, not_function = TRUE)
  } else {
    solver <- find_ilp_solver()
  }
  ilp <- anticlustering_ilp(x, K = 0, FALSE) # k is irrelevant
  solution <- solve_ilp_diversity(ilp, "max", solver)
  ilp_to_groups(solution, nrow(x))
}

# param X is for repetitions
heuristic_wce <- function(X, x) {
  N <- nrow(x)
  rnd_order <- sample(N)
  results <- .C(
    "wce_heuristic", 
    as.double(x[rnd_order, rnd_order]),
    as.integer(N),
    clusters = as.integer(0:(N-1)),
    mem_error = as.integer(0),
    PACKAGE = "anticlust"
  )
  if (results[["mem_error"]] == 1) {
    stop("Could not allocate enough memory.")
  }
  clusters <- results[["clusters"]] + 1
  order_cluster_vector(clusters[order(rnd_order)])
}

#' Compute similarity values for the clique partitioning problem
#' 
#' @param x A data frame / matrix of categorical (nomial / ordinal) attributes
#' @param NAs_equal_each_other If two persons have missing values on the same
#' attribute, does this count as agreement? Defaults to FALSE, in which case
#' these attributes are ignored when computing the sum of agreements and disagreements.
#' 
#' @return An object of class \code{dist} representing pairwise similarity (not (DIS)similartiy!)
#' 
#' @details
#' "the [...] values represent the number of attributes on which vertices
#' [...] disagree, minus the number of attributes on which they agree." 
#' (Brusco et al. 2009, p. 689) By default, NAs are not thought to be equal
#' to each other, as in \code{NA == NA}. Adjust the argument \code{NAs_equal_each_other}
#' if two missing values indicate similarity.
#' 
#' @references 
#' 
#' Brusco, M. J., & Köhn, H. F. (2009). Clustering qualitative data
#' based on binary equivalence relations: neighborhood search
#' heuristics for the clique partitioning problem. Psychometrika, 74,
#' 685–703.
#' 
#' Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm
#' for a clustering problem. Mathematical Programming, 45, 59-96.
#' 
#' @export
cpp_similarities <- function(x, NAs_equal_each_other = FALSE) {
  N <- nrow(x)
  M <- ncol(x)
  output <- rep(0, choose(N, 2))
  x <- apply(x, 2, to_numeric) # by column, convert to integer
  if (NAs_equal_each_other) {
    x[is.na(x)] <- max(x, na.rm = TRUE) + 1 # here: NA is distinct category
  } else {
    x[is.na(x)] <- -1 # In the C implementation, NAs are recognized as -1
  }
  results <- .C(
    "cpp_similarities_", 
    as.integer(x),
    as.integer(N),
    as.integer(M),
    output = as.integer(output),
    PACKAGE = "anticlust"
  )
  lower_tri_to_dist(results$output, N)
}

# convert a vector, filling only the lower triangular of a matrix to dist
# the vector is filled by column (column major I guess)
lower_tri_to_dist <- function(x, N) {
  ret <- matrix(ncol = N, nrow = N)
  ret[lower.tri(ret)] <- x
  as.dist(ret)
}

