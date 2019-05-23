
#' Heuristic anticlustering using random search methods
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used.
#' @param K How many anticlusters should be created.
#' @param preclusters An optional vector representing the preclustering
#'     of the elements in \code{features}
#' @param objective The objective to be maximized, either "distance" or
#'     "variance".
#' @param nrep The number of repetitions tried when assigning elements
#'     to anticlusters when the method is "sampling" or "annealing".
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not used. A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}.
#' @param categories A vector, data.frame or matrix that represents
#'     one or several categorical constraints.
#'
#' @return A vector representing the anticlustering.
#'
#' @noRd
#'

heuristic_anticlustering <- function(features, K, preclusters, objective,
                                     nrep, distances, categories) {

  ## What was the input: features or distances
  use_distances <- FALSE
  if (argument_exists(features)) {
    input <- features
  } else {
    input <- distances
    use_distances <- TRUE
  }

  ## Determine plan for random sampling
  sampling_plan <- "unrestricted"
  if (argument_exists(preclusters)) {
    sampling_plan <- "preclustering"
  } else if (argument_exists(categories)) {
    sampling_plan <- "categorical"
    categories <- merge_into_one_variable(categories)
  }

  dat <- sort_by_group(input, preclusters, categories)
  best_assignment <- random_sampling(dat, K, objective, nrep,
                                     sampling_plan, use_distances)

  ## Return anticluster assignment in original order
  dat[, 1] <- best_assignment
  dat <- sort_by_col(dat, 2)
  anticlusters <- dat[, 1]
  names(anticlusters) <- NULL
  anticlusters
}

#' Merge several grouping variable into one
#'
#' @param categories A vector, data.frame or matrix that represents
#'     one or several categorical constraints.
#'
#' @return A vector representing the group membership (or the combination
#'     of group memberships) as one variable
#'
#' @noRd
#'

merge_into_one_variable <- function(categories) {
  categories <- data.frame(categories)
  factor(do.call(paste0, as.list(categories)))
}

#' Sort data by a grouping variable
#'
#' @param input A data matrix representing features or distances.
#' @param preclusters A vector of precluster affiliations. Can be NULL,
#'     see Details.
#' @param categories A vector that represents categorical constraints.
#'     Can be NULL, see Details.
#'
#' @return An extended data matrix. The first column indicates the
#'     group category of each element (precluster or categorical variable).
#'     The matrix is sorted by the group, i.e. by the first column.
#'     The second column is a unique number identifying the original
#'     position of each element. This column is used in
#'     \code{heuristic_anticlustering} to restore the original order of
#'     the data. The other columns represent the original data inpute that
#'     was passed to this function.
#'
#'  @details
#'
#'  This function sorts the input table by precluster affiliation
#'  or by a categorical variable; neither needs to be present, if no
#'  grouping restrictions are passed (precluster or categories), the data
#'  is not sorted.
#'
#' @noRd
#'

sort_by_group <- function(input, preclusters, categories) {
  sort_by <- 1 # default: data does not need to be sorted
  if (argument_exists(preclusters)) {
    sort_by <- preclusters
  } else if (argument_exists(categories)) {
    sort_by <- categories
  }
  dat <- cbind(sort_by, 1:nrow(input), input)
  dat <- sort_by_col(dat, 1)
  return(dat)
}
