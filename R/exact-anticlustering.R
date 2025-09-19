
#' Solve balanced K-anticluster editing exactly using ILP
#'
#' @param distances A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. Cannot an be an
#'     object of class \code{dist} (i.e., has been converted to matrix
#'     before this function is called)
#' @param K How many anticlusters should be created.
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created.
#'
#' @return A vector representing the anticluster affiliation of
#'     elements.
#' @importFrom utils combn
#'
#' @noRd

exact_anticlustering <- function(data, K, preclustering, cannot_link) {
  
  distances <- convert_to_distances(data)
  N <- nrow(distances)

  if (preclustering == TRUE) {
    ilp <- anticlustering_ilp(distances, N / K)
    solution <- solve_ilp(ilp, "min")
    preclusters <- ilp_to_groups(solution, N)
    ## Fix distances - ensures that the most similar items are assigned
    ## to different groups
    distances <- edit_distances(distances, preclusters, value = (sum(distances) + 1) * (-1))
    ## Edit ILP - objective function and group sizes
    ilp$obj_function <- vectorize_weights(distances)$costs
    ilp$rhs <- c(rep(1, choose(N, 3) * 3),
                 rep((N / K) - 1, N))
    ## Solve edited ILP
    solution <- solve_ilp(ilp)
    assignment <- ilp_to_groups(solution, N)
    return(assignment)
  }
  
  if (argument_exists(cannot_link)) {
    distances[cleanup_cannot_link_indices(cannot_link)] <- -(sum(distances) + 1)
  }

  ## Here the ILP is created without adjusting distances; i.e., true
  ## exact anticlustering
  ilp <- anticlustering_ilp(distances, K)
  solution <- solve_ilp(ilp)
  ilp_to_groups(solution, N)
}

# Ensure that a distance matrix is passed
convert_to_distances <- function(data, squared = FALSE) {
  if (!is_distance_matrix(data)) {
    distances <- dist(data)
    if (squared) distances <- distances^2
  } else {
    distances <- data
  }
  as.matrix(distances)
}
