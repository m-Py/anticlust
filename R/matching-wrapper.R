
#' Matching
#'
#' Conduct K-partite and unrestricted minimum distance matching to find pairs 
#' or larger groups of similar elements.
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
#' @param p The size of the groups; the default is 2, in which case the
#'     function returns pairs.
#' @param groups An optional vector inducing grouping restrictions. If the
#'     vector is passed, the argument \code{p} is ignored and matches are
#'     sought between elements of different groups.
#'
#' @return An integer vector encoding the matches.
#'
#'
#' @details
#' 
#' If the argument \code{features} is passed, matching is done based on the Euclidean
#' distance between data points. To find matches, this function uses
#' the same algorithm as implemented in \code{\link{balanced_clustering}}
#' but it differs with regard to the interface: This function specifies the size of the groups,
#' \code{\link{balanced_clustering}} specifies 
#' the number of the clusters. Moreover, this function makes it possible to specify grouping 
#' restrictions using the \code{groups} argument, thus enabling K-partite minimum matching.
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @examples
#'
#' # Find triplets
#' N <- 300
#' lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
#' triplets <- matching(lds, p = 3)
#' plot_clusters(lds, clustering = triplets, within_connection = TRUE)
#'
#'
#' # Match between different plant species
#' species <- iris$Species != "setosa"
#' matched <- matching(iris[species, 1], groups = iris[species, 5])
#' plot_clusters(
#'   data.frame(Species = as.numeric(iris[species, 5]), Sepal.Length = iris[species, 1]),
#'   clustering = matched,
#'   within_connection = TRUE
#' )
#' 
#' @export
#'

matching <- function(features = NULL, distances = NULL, p = NULL, groups = NULL) {
  data <- process_input(features, distances)
  N <- nrow(data)
  # augment data if input is imbalanced
  augmented <- augment_data(data, p, groups)
  data <- augmented$data
  groups <- augmented$groups
  # pass dummy argument that informs the called functions of the presence
  # of data augmentation
  dummy <- rep(c(FALSE, TRUE), c(N, nrow(data) - N))
  cl <- nn_centroid_clustering(data, p, groups, dummy = dummy)
  cl[1:N] # remove augmented data points
}

# If the data is imbalanced: Augment it
augment_data <- function(data, p, groups) {
  if (argument_exists(groups)) {
    groups <- to_numeric(groups)
    return(augment_kpartite(data, groups))
  } else {
    return(list(data = augment_unrestricted(data, p), groups = NULL))
  }
}

augment_kpartite <- function(data, groups) {
  tab <- table(groups)
  # case 1: groups are balanced; then simply return
  if (all(tab == tab[1])) {
    return(list(data = data, groups = groups))
  }
  # case 2: data must be augmented
  left_over <- max(tab) - tab
  data <- add_dimensions(data, sum(left_over))
  # augment grouping vector
  groups <- c(groups, rep(1:max(groups), left_over))
  list(data = data, groups = groups)
}

# Case: unrestricted matching, data cannot be split in groups of size p
augment_unrestricted <- function(data, p) {
  # how much do we need to add:
  offset <- nrow(data) %% p
  # invert offset (we need to *add* columns)
  offset <- ifelse(offset == 0, 0, offset * (-1) + p)
  data <- add_dimensions(data, offset)
  data
}

# param n: how many data points should be added to the matrix
add_dimensions <- function(data, n) {
  if (is_distance_matrix(data)) {
    # add rows and columns in distance matrix
    data <- add_dimensions_(data, nrow(data) + n, nrow(data) + n)
  } else {
    # only add rows for feature matrix
    data <- add_dimensions_(data, nrow(data) + n, ncol(data))
  }
  data
}

# Add rows / columns to matrix
add_dimensions_ <- function(data, rows, cols) {
  old_dims <- dim(data)
  # add high value -- added values must be dissimilar to all other values!
  extreme_value <- round(max(data) * 100000)
  new_mat <- matrix(extreme_value, nrow = rows, ncol = cols)
  ## add data into augmented matrix
  new_mat[1:old_dims[1], 1:old_dims[2]] <- data
  # if a feature matrix is passed, make half of the new values extreme into the other direction
  if (!is_distance_matrix(data) && rows > old_dims[1]) {
    start <- old_dims[1] + 1
    end <- old_dims[1] + floor((rows - old_dims[1]) / 2)
    new_mat[start:end, ] <- extreme_value * (-1)
  }
  new_mat
}
