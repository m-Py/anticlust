## function centroid_anticlustering()
#
# computes the centroid of all available items and iterates over either
# the items farthest or closest to it, to assign it with the k-1
# items next to the desired number of categories.
#
# algorithm suggested and implemented by m.eik michalke <meik.michalke@hhu.de>
#
# data: matrix or data frame to calculate the distance matrix from
# k: number of target categories
# as_vector: if TRUE, returns a vector with category numbers for each item,
#   otherwise a matrix mit k columns listing the item numbers per category
# forward: if TRUE, always starts with the item closest to the centroid,
#   otherwise the farthest is used as the starting point
#
# examples:
# ac_data <- iris[,1:2]
# x_v <- centroid_anticlustering(ac_data, k=3)
#
# x <- centroid_anticlustering(ac_data, k=3, as_vector=FALSE)
# plot(ac_data[x[,1],], pch=19, col="green")
# points(ac_data[x[,2],], pch=3, col="red")
# points(ac_data[x[,3],], pch=4, col="blue")
# cbind(
#   colMeans(ac_data[x[,1],]),
#   colMeans(ac_data[x[,2],]),
#   colMeans(ac_data[x[,3],])
# )
#
# # reverse the process
# x_fw <- centroid_anticlustering(ac_data, k=3, as_vector=FALSE, forward=TRUE)
# plot(ac_data[x_fw[,1],], pch=19, col="green")
# points(ac_data[x_fw[,2],], pch=3, col="red")
# points(ac_data[x_fw[,3],], pch=4, col="blue")
# cbind(
#   colMeans(ac_data[x_fw[,1],]),
#   colMeans(ac_data[x_fw[,2],]),
#   colMeans(ac_data[x_fw[,3],])
# )
centroid_anticlustering <- function(
  data = NULL,
  distances = NULL,
  k=2,
  as_vector=TRUE,
  forward=FALSE
){
  if (argument_exists(data)) {
    # append column means, coordinates of controid for all dimensions
    data_plus <- rbind(data, colMeans(data))
    # calculate distance matrix including
    # make it a true matrix for easier indexing
  } else if (argument_exists(distances)) {
    # determine a centroid item if the input was a distance matrix
    distances <- as.matrix(distances)
    maxima <- apply(distances, 1, max)
    center_item <- which.min(maxima)
    centroid <- distances[center_item, ]
    dm <- rbind(distances, centroid)
    dm <- cbind(dm, c(centroid, 0))
  }
  cases <- nrow(dm) - 1
  # cases + 1 is now the index number for the centroid variable
  c_idx <- cases + 1
  # get the last row which now holds the distances of all other rows to the centroid
  dist_centroid <- dm[c_idx,]

  # results are stored in a vector; results_k stores the category number for each
  # item, whereas results stores the items ordered by categories
  results <- c()
  results_k <- rep(0, cases)
  names(results_k) <- rownames(dm)[-c_idx]
  ignore <- c(c_idx)
  offset <- 1
  # iteration always starts from the point farthest from/closest to the reference
  # of all points still left. from there, the closest k points are appended to the
  # result vector. this is repeated until no more points are left
  while(length(results) < cases){
    farthest <- next_points(
      m=dist_centroid,
      ignore=c(
        c_idx,
        results
      ),
      max=!forward
    )

    append_to_results <- c(
      farthest,
      next_points(
        m=dm[farthest,],
        ignore=c(
          c_idx,
          farthest,
          results
        ),
        n=k-1,
        max=FALSE
      )
    )

    if(length(results) > k){
      append_to_results <- append_to_results[shift(k=k, start=offset)]
    } else {}

    results <- c(
      results,
      append_to_results
    )

    if(isTRUE(as_vector)){
      results_k[append_to_results] <- 1:k
    } else {}

    if(offset < k){
      offset <- offset +1
    } else {
      offset <- 1
    }
  }

  if(isTRUE(as_vector)){
    return(results_k)
  } else {
    # order the results into a nice matrix
    return(matrix(results, ncol=k, byrow=TRUE))
  }
} ## end function centroid_anticlustering()


## helper function next_points()
# m: a row from the distance matrix (reference point)
# ignore: vector of points (rows/cols) to ignore because already used
#     should at least remove the reference point, otherwise
#     the min criterium returns the reference itself
# n: number of points to return
# min: search for points closest to (min) or farthest from (max) reference?
next_points <- function(
  m,
  ignore,
  n=1,
  max=TRUE
){
  m_sorted <- sort(m[!names(m) %in% ignore], decreasing=max)
  return(names(m_sorted[1:n]))
} ## end function next_points()


## helper function shift()
shift <- function(
  k,
  start
){
  if(start > 1){
    return(c(start:k,1:(start-1)))
  } else {
    return(1:k)
  }
} ## end function shift()
  as.matrix(dist(data_plus))
