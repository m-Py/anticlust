
# Sort a matrix or data frame by one column
#
# @param x A matrix or data frame
# @param col An numeric index or character vector of length 1 to select
#     the column by which the table is sorted.
#
# @return The matrix or data frame sorted by `col`
#
sort_by_col <- function(x, col) {
  if (!class(x) %in% c("data.frame", "matrix"))
    stop("The data that is to be sorted has false structure")
  if (!mode(col) %in% c("character", "numeric"))
    stop("Sorting index has false data type")
  if (length(col) != 1)
    stop("False input, can only select by one column")
  x[order(x[, col]), ]
}

