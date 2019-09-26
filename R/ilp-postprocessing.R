
# Return a clustering vector based on ILP solution for cluster editing
#
# @param solution The solution of the instance returned by
#     `solve_ilp`
#
# @return A clustering vector
#
ilp_to_groups <- function(solution, N) {
  selection_matrix <- matrix(ncol = N, nrow = N)
  selection_matrix[upper.tri(selection_matrix)] <- solution$x
  selection_matrix <- as.matrix(as.dist(t(selection_matrix)))
  diag(selection_matrix) <- 0
  clusters_from_selection_matrix(selection_matrix)
}


# Edit distances
#
# Based on a clustering, distances between elements of the same cluster
# are set to a specified value. This is used to enforce or forbid pairs
# of elements to be part of the same anticluster when the distance
# matrix is passed to the ILP solver.  By considering a preclustering
# of items, very similar items will not be assigned to the same group
# when the fixed distance object is used to create the ILP formulation
# of the item assignment instance.
#
# @param distances A N x N matrix of between-item-distances.
#
# @param clustering A vector representing a clustering; objects that
#     are part of the same cluster are assigned a new distance.
# @param value The value that is assigned to the fixed distances.
#     Defaults to -1,000,000.
#
# @return A distance object containing the fixed distances. For items
#     that are part of the same cluster.
#
#
edit_distances <- function(distances, clustering, value = -1000000) {
  n_groups <- length(unique(clustering))
  for (i in 1:n_groups) {
    items <- which(clustering == i) ## which items are in the i'th cluster
    ## two for-loops to tap all pairwise distances that may require
    ## fixing (a lot of unnecessary iterations, probably)
    for (j in 1:(length(items) - 1)) {
      for (t in 2:length(items)) {
        distances[items[j], items[t]] <- value
        distances[items[t], items[j]] <- value
      }
    }
  }
  return(distances)
}
