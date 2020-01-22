
#' Means and standard deviations by group variable formatted in table
#' 
#' @param features A data frame of features
#' @param groups A grouping vector
#' @param decimals The number of decimals
#' @param na.rm Should NAs be removed prior to computing stats 
#'     (Default = FALSE)
#' @param return_diff Boolean. Should an additional row be printed that 
#'     contains the difference between minimum and maximum
#' 
#' @return A table that illustrates means and standard deviations (in brackets)
#' 
#' @examples 
#'
#' data(iris)
#' mean_sd_tab(iris[, -5], iris[, 5])
#' 
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @export

mean_sd_tab <- function(features, groups, decimals = 2, na.rm  = FALSE, 
                        return_diff = FALSE) {
  K <- length(unique(groups))
  # Means
  means <- by(features, groups, colMeans, na.rm = na.rm)
  rows <- names(means)
  means <- t(matrix(unlist(means), ncol = K))
  mode(means) <- "numeric"
  ncols <- K
  if (return_diff) {
    rows <- c(rows, "Diff")
    ncols <- K + 1
    means <- rbind(means, apply(means, 2, range_diff, na.rm = na.rm))
  }
  means <- apply(means, 2, force_decimals, decimals = decimals)
  # SDS
  sds <- by(features, groups, function(x) apply(x, 2, sd, na.rm = na.rm))
  sds <- t(matrix(unlist(sds), ncol = K))
  mode(sds) <- "numeric"
  if (return_diff) {
    sds <- rbind(sds, apply(sds, 2, range_diff, na.rm = na.rm))
  }
  sds <- apply(sds, 2, force_decimals, decimals = decimals)
  # Merge Means and SDS
  values <- paste0(paste(means, sds, sep = " ("), ")")
  dim(values) <- c(ncols, length(means) / (ncols))
  colnames(values) <- colnames(features)
  rownames(values) <- rows
  values
}

# Determine the maximum range in a vector
range_diff <- function(x, na.rm) {
  x <- round(x, 2) # so that it is consistent in the table
  diff(range(x, na.rm = na.rm))
}

#' Force printing a specified number of decimals for a number
#'
#' @param x the numeric values to be printed
#' @param decimals how many decimals are to be printed. Defaults to 2.
#'
#' @return The number in the required format
#'
#' @examples
#' 
#' force_decimals(c(1.23456, 0.873, 2.3456))
#' force_decimals(c(1.23456, 0.873, 2.3456), 3)
#' 
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @noRd
force_decimals <- function(x, decimals = 2) {
  format(round(x, decimals), nsmall = decimals, scientific = FALSE)
}
