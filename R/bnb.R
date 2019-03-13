
## Branch and bound for anticlustering. Maximizes the within-clusters distances.

## Some global variables that are necessary
# - The problem tree (a list as a first try) `tree`
# - The distances `distances`
# - The current maximum minimum `minmax`
# - n_problems, the length of tree
# - the best anticluster affiliations
# - K the number of anticlusters to be created
# - N the total number of elements

N <- 8
distances <- as.matrix(dist(rnorm(N)))
minmax <- -Inf
n_problems <- 0 # A counter of problems
tree <- list()
## Add first subproblem; first element must be in anticluster 1 (Why? -> See my work on partition enumeration;
## if the first anticluster were > 1, I would obtain redundant partitions)
append_subproblem(NULL, 1, 0)
anticlusters <- NULL
K <- 2


#' An element of the tree is a problem
#'
#' @param anticlusters The anticluster affiliations of the old elements
#' @param anticluster The anticluster affiliation of the new element
#' @param objective The objective based on the old elements
#'
#' A problem (-> a list) is appended to `tree`. There is no return value.
#'
append_subproblem <- function(anticlusters, new, objective) {
  problem <- list(clusters = c(anticlusters, new), prev_value = objective)
  tree[[n_problems + 1]] <<- problem
  n_problems <<- n_problems + 1
}

#' Compute objective for the current problem
#'
#' Always works on `tree[[1]]`. Uses the "previous" objective
#' plus the objective that is given
#'
#'
compute_objective <- function() {
  ## Select all elements within the same cluster as the newly added element
  anticlusters <- tree[[1]]$clusters
  current_element <- length(anticlusters)
  current_cluster <- anticlusters[current_element]
  other_elements <- setdiff(which(anticlusters == current_cluster), current_element)
  ## Get all connections of current element.
  connections <- expand.grid(other_elements, current_element)
  distances_of_connections <- distances[as.matrix(connections)]
  ## Add sum of new distances to previous objective
  tree[[1]]$prev_value <<- tree[[1]]$prev_value + sum(distances_of_connections)
}




#' Process a subproblem:
#'
#' - compute objective
#' - decide if objective still good enough (TODO)
#' - add split into new subproblems if appropriate
process_problem <- function() {
  objective <- compute_objective()
  ## replace current minmax by objective if it is better
  if (objective > minmax) {
    minmax <<- objective
  }
  ## Compute the maximum that can still be achieved:

  ## - Pretend that all of the to follow elements are in the same cluster (this is a
  ## very bad heuristic, make better)
  remaining_elements <- (length(tree[[1]]$clusters) + 1):N
  ## distances between all remaining elements
  all_remaining <- distances[remaining_elements,][lower.tri(distances[remaining_elements,])]
  best_possible <- sum(all_remaining)

  ## is best possible value worse than the minimum maximum to be expected?
  if (best_possible >= minmax) {
    for (k in 1:K) {
      ## check that each cluster occurs legally often
      if (sum(tree[[1]]$clusters == k) <= N / K) {
        append_subproblem(tree[[1]]$clusters, k, objective)
      }
    }
  }
  tree[[1]] <<- NULL
  n_problems <<- n_problems - 1
}


## Get all together
for (i in 1:4) {
  cat(" \n\n ---------- \n\n")
  print(tree)
  process_problem()
}

