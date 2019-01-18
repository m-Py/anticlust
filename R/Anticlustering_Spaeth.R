
## Try implementing Spaeth's (1987) algorithm

# Uses an "exchange method":

# "That method, applied for (2), improves some (random) initial
# partition by successively moving on trial each object from
# its cluster to all the other ones, and by shifting it, if there is
# any reduction at all, to that one where the first term on the right
# side of (4) decreases most, otherwise taking the next object and
# finally passing through all the objects until no further improvement
# occurs.


## all is chaos:

N <- 10
items <- as.matrix(data.frame(rnorm(N), rnorm(N)))
#distances <- as.matrix(dist(items))
anticlusters <- rep(1:2, N/2)

squared_distances <- dist_from_centers(items, matrix(means, ncol = 1))
obj_value(anticlusters, squared_distances)

as.matrix(compute_centers(items, anticlusters))

#' Get feature centers given anticluster grouping
compute_centers <- function(items, anticlusters) {
  by(as.matrix(items), anticlusters, colMeans)
}

obj_value <- function(anticlusters, squared_distances) {
  ## determine distances within each group
  summed_distances <- by(squared_distances, INDICES = anticlusters, sum)
  ## determine objective as the sum of all distances per group
  return(sum(summed_distances))
}

#' Check if a data point should be shifted
#'
#' @return The anticluster the data point should be shifted to
shift <- function(point, items, anticlusters, current_best) {
  items <- matrix(items)
  ## Select all anticlusters the current point is not part of
  other_clusters <- setdiff(anticlusters, point)
  for (i in seq_along(other_clusters)) {
    # flip item and check if objective is improved
    tmp <- items
    tmp
  }
}

shift(1, items, anticlusters, obj_value(anticlusters, squared_distances))
