
#' Merge several grouping variable into one
#'
#' @param categories A vector, data.frame or matrix that represents
#'     one or several categorical constraints.
#'
#' @return A vector representing the group membership (or the combination
#'     of group memberships) as one variable
#'
#' @noRd
#'

merge_into_one_variable <- function(categories) {
  if (is.null(categories)) {
    return(NULL)
  }
  categories <- data.frame(categories)
  categories <- factor(do.call(paste0, as.list(categories)))
  # sort as numeric to get consistent return value
  order_cluster_vector(to_numeric(categories))
}

#' Converting any categorical input to integer
#' 
#' @param x a vector
#' 
#' @return A vector of \code{length(x)} where each element is integer
#' 
#' @examples 
#' iris$Species
#' to_numeric(iris$Species)
#' 
#' @noRd
#' 
to_numeric <- function(x) {
  as.numeric(factor(x))
}

# Order a clustering vector
# 
# For a clustering vector, ensure that the first cluster that occurs
# in the vector is 1, the next 2, etc ... until K
#
# param clusters: A clustering vector with elements 1, ..., K indicating cluster membership
# return: A clustering vector in order
order_cluster_vector <- function(clusters) {
  unique_clusters <- unique(clusters)
  # deal with NA
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]
  K <- length(unique_clusters)
  clusters <- factor(clusters, levels = unique_clusters, labels = 1:K)
  as.numeric(clusters)
}
