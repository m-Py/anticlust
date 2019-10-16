
#' Solve balanced K-anticluster editing exactly using ILP
#'
#' @param distances A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. Cannot an be an
#'     object of class \code{dist} (i.e., has been converted to matrix
#'     before this function is called)
#' @param K How many anticlusters should be created.
#' @param solver A string identifying the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created.
#'
#' @return A vector representing the anticluster affiliation of
#'     elements.
#'
#' @noRd

exact_anticlustering <- function(distances, K, solver, preclustering) {

  # compute distance matrix if it was not passed
  if ("features" %in% class(distances)) {
    distances <- as.matrix(dist(distances))
  }

  N <- nrow(distances)

  if (preclustering == TRUE) {
    ilp <- anticlustering_ilp(distances, N / K,
                              solver = solver)
    solution <- solve_ilp(ilp, solver, "min")
    preclusters <- ilp_to_groups(solution, N)
    ## Fix distances - ensures that the most similar items are assigned
    ## to different groups
    distances <- edit_distances(distances, preclusters, value = (sum(distances) + 1) * (-1))
    ## Edit ILP - objective function and group sizes
    ilp$obj_function <- vectorize_weights(distances)$costs
    ilp$rhs <- c(rep(1, choose(N, 3) * 3),
                 rep((N / K) - 1, N))
    ## Solve edited ILP
    solution <- solve_ilp(ilp, solver)
    assignment <- ilp_to_groups(solution, N)
    return(assignment)
  }

  ## Here the ILP is created without adjusting distances; i.e., true
  ## exact anticlustering
  ilp <- anticlustering_ilp(distances, K, solver = solver)
  solution <- solve_ilp(ilp, solver)
  ilp_to_groups(solution, N)
}
