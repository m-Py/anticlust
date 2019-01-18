
## Work on implementation of anticlustering

# First: implement hierarchical procedure by Valev (1998)


# 1.  Each element is its own anticluster
# 2. For each anticluster, find the minimum distance between one of its
#    elements and the element of a different anticluster
# 3. Merge the anticlusters where the index from 2 is maximum
# 4. Repeat 2 & 3 until the required number of anticlusters is found

anti <- data.frame(value = rnorm(10), anticluster = 1:10)
distances <- dist(anti$value)


#' Minimum distances to elements of other anticlusters
#'
#' @param distances A matrix or dist object representing distances
#'   between all elements
#' @param anticlusters A numeric vector indicating anticluster membership
#'   for each element
#'
#' @return A data.frame containing information for each element on the
#'   closest element in a different anticluster
#'
#' @example
#' distances <- dist(rnorm(10))
#' mindist_elementwise(distances, rep(1:5, 2))
mindist_elementwise <- function(distances, anticlusters) {
  ## Some error handling
  if (nrow(distances) != length(anticluster))
    stop("Error: Distances and anticlusters indicate a different number of elements.")
  distances <- as.matrix(distances)
  N <- nrow(distances)
  # Initialize vector to store minima for each anticluster
  minima <- matrix(nrow = N, ncol = 4)
  colnames(minima) <- c("minimum", "anticluster", "target_anticluster", "target_element")
  for (i in 1:N) {
    # 1. Find out which anticluster the i'th element is in
    anticluster <- anticlusters[i]
    # 2. Select distances to other elements for element i
    distance <- distances[i, ]
    # 3. Determine which elements are in different anticlusters than i
    in_different_anticluster <- anticlusters != anticluster
    # 4. Get minimum distance for i'th element
    minimum <- min(distance[in_different_anticluster])[1] # only select one target!
    # 5. Get target: Which element is closest to i'th element?
    target_element <- which(distance == minimum & in_different_anticluster)[1]
    ## Assert that a target element was found
    if (length(target_element) != 1 | is.na(target_element))
      stop("Something went wrong: Target element was not found")
    # 6. Fill values into matrix
    minima[i, ] <- c(minimum, anticluster, anticlusters[target_element], target_element)
  }
  minima <- data.frame(minima)
  ## Assert that minimum was found for a different anticluster
  if (any(minima$anticluster == minima$target_anticluster))
    stop("Something went wrong: Target anticluster cannot be the same as the original anticluster")
  return(minima)
}

#' Minimum distances per anticluster
#'
#' @param distances A matrix or dist object representing distances
#'   between all elements
#' @param anticlusters A numeric vector indicating anticluster membership
#'   for each element
#'
#' @return Information on the anticlusters to be merged
mindist_anticluster <- function(distances, anticlusters) {
  ## Compute minimum distances for each element
  minima <- mindist_elementwise(distances, anticlusters)
  ## Find minimum per anticluster
  minima <- by(minima, minima$anticluster, function(x) x[x$minimum == min(x$minimum), ])
  ## Ensure that per anticluster, only only one row is returned:
  minima <- do.call(rbind, lapply(minima, function(x) x[1, ]))
  if (nrow(minima) != length(unique(anticlusters)))
    stop("Somethin went wrong, there is not a minimum for each anticluster")
  return(minima)
}


#' Information which clusters should be merged
#'
#' @param distances A matrix or dist object representing distances
#'   between all elements
#' @param anticlusters A numeric vector indicating anticluster membership
#'   for each element
#'
#' @return A vector of length 2 indicating the anticlusters to be merged.
merge_which <- function(distances, anticlusters) {
  minima <- mindist_anticluster(distances, anticlusters)
  max_min <- max(minima$minimum)
  where_maxmin <- which(minima$minimum == max_min)[1] # one target
  unlist(minima[where_maxmin, c("anticluster", "target_anticluster")])
}

distances <- dist(rnorm(8))
anticlusters <- c(rep(2:4, 2), 2, 2)
mindist_elementwise(distances, anticlusters)
mindist_anticluster(distances, anticlusters)
merge_which(distances, anticlusters)
