
#' Matching
#'
#' Conduct K-partite or unrestricted (minimum distance) matching to
#' find pairs or groups of similar elements. By default, finding
#' matches is based on the Euclidean distance between data points, but
#' a custom distance measure can also be employed.
#'
#' @param features A numeric vector, matrix or data.frame of data
#'     points.  Rows correspond to elements and columns correspond to
#'     features. A vector represents a single numeric feature.
#' @param distances Alternative data argument that can be used if
#'     \code{features} is not passed. An N x N matrix representing the
#'     pairwise dissimilarities between N elements. Larger values
#'     indicate higher dissimilarity. Can be an object of class
#'     \code{dist} (e.g., returned by \code{\link{dist}} or
#'     \code{\link{as.dist}}) or a \code{matrix} where the entries of
#'     the upper and lower triangular matrix represent the pairwise
#'     dissimilarities.
#' @param p The size of the groups; the default is 2, in which case
#'     the function returns pairs.
#' @param match_between An optional vector, \code{data.frame} or
#'     matrix representing one or several categorical constraints. If
#'     passed, the argument \code{p} is ignored and matches are sought
#'     between elements of different categories.
#' @param match_within An optional vector, \code{data.frame} or matrix
#'     representing one or several categorical constraints. If passed,
#'     matches are sought between elements of the same category.
#' @param match_extreme_first Logical: Determines if matches are first
#'     sought for extreme elements first or for central
#'     elements. Defaults to \code{TRUE}.
#' @param target_group lawl
#'
#' @return An integer vector encoding the matches. See Setails for
#'     more information.
#'
#'
#' @details
#' 
#' If the argument \code{features} is passed, matching is done based
#' on the Euclidean distance between data points. If the argument
#' \code{distances} is passed, entries in this matrix define
#' dissimilarity between data points. To find matches, the algorithm
#' proceeds by selecting a target element and then searching its
#' nearest neighbours. Critical to the behaviour or the algorithm is
#' the order in which target elements are selected. By default the
#' most extreme elements are selected first, i.e., elements with the
#' highest distance to the center of the data set (see
#' \code{\link{balanced_clustering}}). By setting the argument
#' \code{match_extreme_first} to \code{FALSE}, it is possible to
#' enforce that elements close to the center are first selected as
#' targets. If the argument \code{match_between} is passed and the
#' groups specified via this argument are of different size, target
#' elements are always selected from the smallest group (because in
#' this group, all elements can be matched).
#' 
#' The output is an integer vector encoding which elements have been
#' matched. The grouping numbers are sorted by similarity. That is,
#' elements with the grouping number »1« are most similar to each
#' other, followed by 2 etc (groups having the same similarity index
#' are still assigned a different grouping number, though). Similarity
#' is measured as the sum of pairwise (Euclidean) distances within
#' groups (see \code{\link{distance_objective}}).  Some elements of
#' the output may be \code{NA}. This happens if it is not possible to
#' evenly split the item pool evenly into groups of size \code{p} or
#' if the categories described by the argument \code{match_between}
#' are of different size; unmatched items are assigned \code{NA}.
#' 
#' @note It is possible to specify grouping restrictions via 
#' \code{match_between} and \code{match_within} at the same time.
#' 
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' # Find triplets
#' N <- 29 # two elements will not be matched
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' triplets <- matching(lds, p = 3)
#' table(triplets)
#' sum(is.na(triplets))
#' plot_clusters(
#'   lds,
#'   clustering = triplets,
#'   within_connection = TRUE
#' )
#'
#' # Bipartite matching with unequal-sized groups:
#' # Only selects matches for some elements
#' N <- 100
#' data <- matrix(rnorm(N), ncol = 1)
#' groups <- sample(1:2, size = N, replace = TRUE, prob = c(0.8, 0.2))
#' matched <- matching(data[, 1], match_between = groups)
#' plot_clusters(
#'   cbind(groups, data), 
#'   clustering = matched, 
#'   within_connection = TRUE
#' )
#' 
#' # Match objects from the same category only
#' matched <- matching(
#'   schaper2019[, 3:6], 
#'   p = 3, 
#'   match_within = schaper2019$room
#' )
#' table(matched, schaper2019$room)
#' 
#' # Match between different plant species in the »iris« data set
#' species <- iris$Species != "versicolor"
#' matched <- matching(
#'   iris[species, 1], 
#'   match_between = iris[species, 5]
#' )
#' # Adjust `match_extreme_first` argument
#' matched2 <- matching(
#'   iris[species, 1], 
#'   match_between = iris[species, 5],
#'   match_extreme_first = FALSE
#' )
#' # Plot the matching results
#' par(mfrow = c(1, 2))
#' data <- data.frame(
#'   Species = as.numeric(iris[species, 5]),
#'   Sepal.Length = iris[species, 1]
#' )
#' plot_clusters(
#'   data,
#'   clustering = matched,
#'   within_connection = TRUE,
#'   main = "Extreme elements matched first"
#' )
#' plot_clusters(
#'   data,
#'   clustering = matched2,
#'   within_connection = TRUE,
#'   main = "Central elements matched first"
#' )
#' par(mfrow = c(1, 1))
#' 
#' 
#' @export
#'

matching <- function(
  features = NULL, 
  distances = NULL, 
  p = 2, 
  match_between = NULL,
  match_within = NULL,
  match_extreme_first = TRUE,
  target_group = NULL
) {
  data <- process_input(features, distances)
  if (argument_exists(target_group)) {
    id <- which(match_between == target_group)[1]
  }
  match_between <- merge_into_one_variable(match_between)
  if (argument_exists(target_group)) {
    target_group <- match_between[id]
  }
  if (argument_exists(match_within)) {
    cl <- match_within(data, p, match_between, match_within, match_extreme_first, target_group)
  } else {
    cl <- nn_centroid_clustering(data, p, match_between, match_extreme_first, target_group)
  }
  # Before returning: order the group numbers by objective - most similar 
  # matches have lower indices
  sort_by_objective(cl, data)
}

# conduct a matching for each category if `match_within` is passed
match_within <- function(data, p, match_between, match_within, match_extreme_first, target_group) {
  match_within <- merge_into_one_variable(match_within)
  N <- nrow(data)
  cl <- rep(NA, N)
  c <- length(unique(match_within))
  for (i in 1:c) {
    tmp_data <- subset_data_matrix(data, match_within == i)
    cl_tmp <- nn_centroid_clustering(
      tmp_data, 
      p, 
      match_between[match_within == i], 
      match_extreme_first
    )
    # ensure that different cluster numbers are given to different groups
    cl[match_within == i] <- ifelse(is.na(cl_tmp), NA, paste0(cl_tmp, "_", i))
  }
  to_numeric(cl)
}

sort_by_objective <- function(cl, data, N) {
  N <- nrow(data)
  selected <- (1:N)[!is.na(cl)]
  cl_sub <- cl[selected]
  cl_sub <- order_cluster_vector(cl_sub)
  if (is_distance_matrix(data)) {
    data <- data[selected, selected]
    objectives <- sapply(
      1:max(cl_sub), 
      function(x) sum(as.dist(data[cl_sub == x, cl_sub == x]))
    )
  } else {
    data <- data[selected, , drop = FALSE]
    objectives <- sapply(
      1:max(cl_sub), 
      function(x) sum(dist(data[cl_sub == x, ]))
    )
  }
  # recode original matching labels according to objective
  N <- length(cl_sub)
  # sorry for this code, but it works and it was really complicated
  one <- data.frame(match = cl_sub, order_matches = 1:length(cl_sub))
  two <- data.frame(match = order(objectives), objective_id = 1:length(objectives))
  merged <- merge(one, two)
  new <- rep(NA, N)
  new[merged$order_matches] <- merged$objective_id
  cl[!is.na(cl)] <- new
  cl
}
