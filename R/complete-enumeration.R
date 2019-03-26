
#' Complete enumeration approach to anticluster editing
#'
#' Computes the optimum anticluster editing objective but is probably
#' very slow because all possible partitions of anticlusters are
#' enumerated.
#'
#' @param features A vector, matrix or data.frame of data points. Rows
#'     correspond to elements and columns correspond to features. A
#'     vector represents a single feature.
#' @param K How many anticlusters should be created.
#'
#' @return The anticluster affiliation
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
#' results <- enum_anticlustering(features, K = K)

enum_anticlustering <- function(features = NULL, distances = NULL, K) {


  if (!argument_exists(features) && !argument_exists(distances)) {
    stop("One of the arguments 'features' or 'distances' must be given.")
  }

  if (argument_exists(features) && argument_exists(distances)) {
    stop("Only pass one of the arguments 'features' or 'distances'.")
  }

  if (argument_exists(features)) {
    distances <- dist(features)
  }
  distances <- as.matrix(distances)

  ## How many items are to be reassigned:
  N <- nrow(distances)
  ## Initialize a vector that encodes the assignment to groups
  anticlusters  <- sort(rep_len(1:K, N))
  ## Initialize objective
  best_objective <- -Inf
  best_assign <- NULL

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
    cur_obj <- distance_objective(distances, anticlusters)
    ## Better fit was found, save the assignment
    if (cur_obj > best_objective) {
      best_assign <- anticlusters
      best_objective <- cur_obj
    }

    anticlusters <- next_permutation(anticlusters)
  }

  return(best_assign)
}


#' Get the first occurrence of each value in a vector, sorted by values
#'
#' @param anticlusters a vector representing anticluster affiliations
#' @param K The number of anticlusters
#' @param first defaults to 1, and therefore returns the first occurence
#'     of each anticluster. Can be changed to return the second, third
#'     etc. occurence.
#'
#' @return The first occurrence of each value 1:K
#'
#' @noRd
#'
first_occurrences <- function(anticlusters, K, first = 1) {
  first_occurrences <- rep(NA, K)
  for (i in 1:K) {
    first_occurrences[i] <- which(anticlusters == i)[first]
  }
  first_occurrences
}
