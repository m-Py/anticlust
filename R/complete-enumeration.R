
#' Complete enumeration approach to anticlustering
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param K How many anticlusters should be created.
#' @param talk Boolean. If `TRUE`, the function will print its progress.
#' @param objective The objective to be maximized, either "distance" or
#'     "variance".
#'
#' @return A \code{list}. The first element is the best assignment.
#'   The second element is the best objective.
#'
#' @export
#'
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' n_features <- 2
#' N <- 10
#' K <- 2
#' features <- matrix(runif(n_features * N), ncol = n_features)
#' results <- enum_anticlustering(features, K, objective = "distance")

enum_anticlustering <- function(features, K, objective = "distance") {

  ## How many items are to be reassigned:
  N <- nrow(features)
  ## Initialize a vector that encodes the assignment to groups
  anticlusters  <- sort(rep_len(1:K, N))
  ## Initialize objective
  best_objective <- get_objective(features, anticlusters, objective)
  best_assign <- anticlusters

  i <- 1
  repeat {
    ## ending condition: next permutation is starting point
    if (anticlusters[1] != 1) { # only redundant partitions from this point
      break
    }

    ## Only compute objective for unique partitions:
    if (is.unsorted(first_occurrences(anticlusters, K))) {
      anticlusters <- next_permutation(anticlusters)
      next
    }

    ## Check the objective value
    cur_obj <- get_objective(features, anticlusters, objective)
    ## Better fit was found, save the assignment
    if (cur_obj > best_objective) {
      best_assign <- anticlusters
      best_objective <- cur_obj
    }

    i <- i + 1
    anticlusters <- next_permutation(anticlusters)
  }

  return(list(anticlusters = best_assign, objective = best_objective))
}


#' Get the first occurrence of each value in a vector, sorted by values
#'
#' @param anticlusters a vector representing anticluster affiliations
#'
#' @return The first occurrence of each value 1:K
#'
#' @noRd
#'
first_occurrences <- function(anticlusters, K) {
  first_occurrences <- rep(NA, K)
  for (i in 1:K) {
    first_occurrences[i] <- which(anticlusters == i)[1]
  }
  first_occurrences
}
