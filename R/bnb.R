
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
  anticlusters <- anticlustering(features, K, method = "sampling", nrep = 2)
  minmax <- get_objective(features, anticlusters, "distance")
  distances <- as.matrix(dist(features))
  n_problems <- 0 # A counter of problems
  tree <- list()
  N <- nrow(features)

  # An element of the tree is a problem
  #
  # @param new_clusters The anticluster affiliations of all elements in
  #     the candidate solution
  # @param objective The objective based on the old elements
  #
  # A problem (-> a list) is appended to `tree`. There is no return value.
  #
  append_subproblem <- function(new_clusters, objective) {
    problem <- list(clusters = new_clusters, prev_value = objective)
    tree[[n_problems + 1]] <<- problem
    n_problems <<- n_problems + 1
  }

  # Compute objective for the current problem
  #
  #
  # @param problem, will always be `tree[[1]]`
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
  # - decide if objective still good enough (TODO)
  # - add split into new subproblems if appropriate
  process_problem <- function() {
    objective <- compute_objective(tree[[1]])
    length_of_candidate <- length(tree[[1]]$clusters)
    ## replace current minmax by objective if it is better
    if (objective >= minmax) {
      minmax <<- objective
      anticlusters <<- tree[[1]]$clusters
    }

    ## Compute the maximum that can still be achieved if there is
    ## not all elements in the candidate solution.
    ## Idea: for each element that is still to be inserted into a
    ## candidate solution, select the maximum N / K-1 distances to elements.
    if (length_of_candidate < N) {
      remaining_elements <- (length_of_candidate + 1):N
      maxima <- c()
      for (i in remaining_elements) {
        ## all distances for element i
        tmp_dists <- sort(distances[i, ], decreasing = TRUE)
        maxima <- c(maxima, sum(tmp_dists[1:(N / K - 1)]))
      }
      best_possible <- sum(maxima) + objective
      ## is best possible value worse than the minimum maximum to be expected?
      if (best_possible >= minmax) {
        for (k in 1:K) {
          ## 1) Check that each cluster occurs legally often in new candidate
          ## 2) Check that the new candidate is not a redundant partition
          ## (this means that k cannot be larger than the length of the new
          ## candidate!)
          new_clusters <- c(tree[[1]]$clusters, k)
          if (sum(new_clusters == k) <= N / K &&
              k <= length_of_candidate + 1) {
            ## new problem gets as "previous objective" the current objective
            append_subproblem(new_clusters, objective)
          }
        }
      }
    }

    ## remove the problem that was processed
    tree[[1]] <<- NULL
    n_problems <<- n_problems - 1
  }

  ## Add first subproblem; first element must be in anticluster 1
  ## (Why? -> See my work on partition enumeration;
  ## if the first anticluster were > 1, I would obtain redundant partitions)

  append_subproblem(1, 0)
  ## Get all together
  while (length(tree) != 0) {
    process_problem()
  }

  return(list(objective = minmax, anticlusters = anticlusters))
}

