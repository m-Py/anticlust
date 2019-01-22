
#' Utility function to sort a matrix or data frame by one column
#'
#' @param x A matrix or data frame
#' @param row An index for row selection (numeric; character; logical)
#'
#' @return The matrix or data frame sorted by `row`
#' @examples
#' mat <- cbind(10:1, 1:10)
#' mat
#' sort_by_col(mat, 1)
sort_by_col <- function(x, col) {
  if (!class(x) %in% c("data.frame", "matrix"))
    stop("The data that is to be sorted has false structure")
  if (!mode(col) %in% c("character", "numeric", "logical"))
    stop("Sorting index has false data type")
  x[order(x[, col]), ]
}

