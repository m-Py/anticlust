
#' Sort a matrix or data frame by one column
#'
#' @param x A matrix or data frame
#' @param row An numeric index or character vector for selecting
#'   the column by which the table is sorted. Needs to have length 1.
#'
#' @return The matrix or data frame sorted by `row`
#'
#' @export
#'
#' @examples
#' mat <- cbind(10:1, 1:10)
#' mat
#' sort_by_col(mat, 1)
sort_by_col <- function(x, col) {
  if (!class(x) %in% c("data.frame", "matrix"))
    stop("The data that is to be sorted has false structure")
  if (!mode(col) %in% c("character", "numeric"))
    stop("Sorting index has false data type")
  if (length(col) != 1)
    stop("False input, can only select by one column")
  x[order(x[, col]), ]
}

