
#' Means and standard deviations by group variable formatted in table
#' 
#' @param features A data frame of features
#' @param groups A grouping vector
#' @param decimals The number of decimals
#' 
#' @return A table that illustrates means and standard deviations (in brackets)
#' 
#' @examples 
#'
#' data(iris)
#' mean_sd_tab(iris[, -5], iris[, 5])
#' 
#' @author: Martin Papenberg
#' @export

mean_sd_tab <- function(features, groups, decimals = 2, na.rm  = FALSE) {
  K <- length(unique(groups))
  # Means
  means <- by(features, groups, colMeans, na.rm = na.rm)
  means <- t(matrix(unlist(means), ncol = K))
  mode(means) <- "numeric"
  means <- rbind(means, apply(means, 2, range_diff, na.rm = na.rm))
  means <- apply(means, 2, prmisc::force_decimals, decimals = decimals)
  # SDS
  sds <- by(features, groups, function(x) apply(x, 2, sd, na.rm = na.rm))
  sds <- t(matrix(unlist(sds), ncol = K))
  mode(sds) <- "numeric"
  sds <- rbind(sds, apply(sds, 2, range_diff, na.rm = na.rm))
  sds <- apply(sds, 2, prmisc::force_decimals, decimals = decimals)
  # Merge Means and SDS
  values <- paste0(paste(means, sds, sep = " ("), ")")
  dim(values) <- c(K + 1, length(means) / (K + 1))
  colnames(values) <- colnames(features)
  values
}

# Determine the maximum range in a vector
range_diff <- function(x, na.rm) {
  x <- round(x, 2) # so that it is consistent in the table
  diff(range(x, na.rm = na.rm))
}
