#' Dispersion
#' 
#' Compute the dispersion objective for a given clustering (i.e., the minimum distance 
#' between two elements within the same cluster)
#' 
#' @noRd

dispersion_objective <- function(clusters, distances) {
  objectives <- sapply(
    sort(unique(clusters)), 
    function(x) min(as.dist(distances[clusters == x, clusters == x]))
  )
  min(objectives)
}
