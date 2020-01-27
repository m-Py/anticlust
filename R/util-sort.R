
#' Sort a matrix or data frame by one column
#'
#' @param x A matrix or data frame
#' @param col An numeric index or character vector of length 1 to select
#'     the column by which the table is sorted.
#' @param decreasing logical. Should the sort order be increasing or
#'     decreasing?
#'
#' @return The matrix or data frame sorted by `col`
#'
#' @noRd
sort_by_col <- function(x, col, decreasing = FALSE) {
  if (length(col) != 1)
    stop("False input, can only select by one column")
  x[order(x[, col], decreasing = decreasing), ]
}
