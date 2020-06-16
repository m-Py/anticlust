
#' Plot similarity objective by cluster
#' 
#' @param x The data input. Can be one of two structures: (1) A data matrix
#'     where rows correspond to elements and columns correspond to
#'     features (a single numeric feature can be passed as a vector). (2)
#'     An N x N matrix dissimilarity matrix; can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities. 
#' @param groups A grouping vector of length N, usually the output
#'     of \code{\link{matching}}.
#' 
#' @return The diversity (sum of distances) by group.
#' 
#' @details Plots the sum of pairwise distances by group.
#' 
#' @export
#' 
#' @examples
#' 
#' # Match elements and plot similarity by match
#' N <- 100
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' pairs <- matching(lds, p = 2)
#' plot_similarity(lds, pairs)
#' 
#' @seealso
#'
#' \code{\link{diversity_objective}}
#'
#' @author 
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 

plot_similarity <- function(x, groups) {
  x <- as.matrix(x)
  select <- !is.na(groups)
  x <- subset_data_matrix(x, select)
  groups <- groups[select]
  diversity <- diversity_objective_by_group(groups, x)
  plot(diversity, ylab = "Sum of distances", xlab = "Group")
  return(invisible(diversity))
}
