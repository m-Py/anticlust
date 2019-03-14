
#' Simple branch and bound approach to anticlustering
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param K How many anticlusters should be created.
#'
#' @return A list containing the objective and the anticluster
#'     affiliations.
#'
#' @export
#'


bnb_anticlustering <- function(features, K) {

  ## Some "global" variables are necessary
  # - tree - The problem tree (a list as a first try)
  # - `distances` - the distances between all elements
  # - `minmax` - The current minimum maximum
  # - n_problems- the length of `tree`
  # - anticluster - the current best anticluster affiliations
  # - K the number of anticlusters to be created
  # - N the total number of elements
  # - temporily: candidates - a list containing all generated candidates
  # - anticlusters - the current best partitioning

  features <- as.matrix(features)
  anticlusters <- anticlustering(features, K, method = "sampling", nrep = 50)
  minmax <- get_objective(features, anticlusters, "distance")
  distances <- as.matrix(dist(features))
  tree <- new.queue()
  N <- nrow(features)

  # Compute objective for the current problem
  #
  #
  # @param problem, will always be `dequeue(tree)`
  #
  # @return The objective value given all elements in the candidate solution
  #
  #
  compute_objective <- function(problem) {
    ## Select all elements within the same cluster as the newly added element
    ## (most recent element, i.e., element at last position in `problem$clusters`
    anticlusters <- problem$clusters
    current_element <- length(anticlusters)
    current_cluster <- anticlusters[current_element]
    other_elements <- setdiff(which(anticlusters == current_cluster), current_element)
    ## Get all connections of current element.
    connections <- expand.grid(other_elements, current_element)
    ## Get all distances between current element and elements of the same anticluster
    distances_of_connections <- distances[as.matrix(connections)]
    ## Add sum of new distances to previous objective
    objective <- problem$prev_value + sum(distances_of_connections)
    return(objective)
  }

  # Process a subproblem:
  #
  # - compute objective
  # - decide if objective still good enough
  # - add split into new subproblems if appropriate

  process_problem <- function() {
    problem <- dequeue(tree) # takes uppermost candidate
    objective <- compute_objective(problem)
    length_of_candidate <- length(problem$clusters)
    ## replace current minmax by objective if it is better
    if (objective >= minmax) {
      minmax <<- objective
      anticlusters <<- problem$clusters
    }

    ## Compute the maximum that can still be achieved if there is
    ## not all elements in the candidate solution.
    ## Idea: for each element that is still to be inserted into a
    ## candidate solution, select the maximum N / K-1 distances to elements.
    if (length_of_candidate < N) {
      remaining_elements <- (length_of_candidate + 1):N
      maxima <- rep(NA, length(remaining_elements))
      for (i in seq_along(remaining_elements)) {
        ## all distances for element i
        tmp_dists <- sort(distances[remaining_elements[i], ], decreasing = TRUE)
        maxima[i] <- sum(tmp_dists[1:(N / K - 1)])
      }
      best_possible <- sum(maxima) + objective
      ## is best possible value worse than the minimum maximum to be expected?
      if (best_possible >= minmax) {
        for (k in 1:K) {
          ## 1) Check that each cluster occurs legally often in new candidate
          ## 2) Check that the new candidate is not a redundant partition
          ## (this means that k cannot be larger than the length of the new
          ## candidate!)
          ## Only append new subproblem if both conditions are satisfied
          new_clusters <- c(problem$clusters, k)
          if (sum(new_clusters == k) <= N / K &&
              k <= length_of_candidate + 1) {
            ## new problem gets as "previous objective" the current objective
            enqueue(tree, list(clusters = new_clusters, prev_value = objective))
          }
        }
      }
    }
  }

  ## Add first subproblem; first element must be in anticluster 1
  ## (Why? -> See my work on partition enumeration;
  ## if the first anticluster were > 1, I would obtain redundant partitions)
  enqueue(tree, list(clusters = 1, prev_value = 0))
  ## Get all together
  while (!is.empty(tree)) {
    process_problem()
  }

  return(list(objective = minmax, anticlusters = anticlusters))
}


## Define a queue in R using environments
## (source: https://www.researchgate.net/post/What_is_the_queue_data_structure_in_R)

new.queue <- function() {
  ret <- new.env()
  ret$front <- new.env()
  ret$front$q <- NULL
  ret$front$prev <- NULL
  ret$last <- ret$front
  return(ret)
}
## add to end of queue
enqueue <- function(queue, add){
  queue$last$q <- new.env()
  queue$last$q$prev <- queue$last
  queue$last <- queue$last$q
  queue$last$val <- add
  queue$last$q <- NULL
}
## return front of queue and remove it
dequeue <- function(queue){
  if (is.empty(queue)) {
    stop("Attempting to take element from empty queue")
  }
  value <- queue$front$q$val
  queue$front <- queue$front$q
  queue$front$q$prev <- NULL
  return(value)
}
is.empty <- function(queue){
  return(is.null(queue$front$q))
}
