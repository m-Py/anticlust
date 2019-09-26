
#' Solve balanced K-anticluster editing exactly using ILP
#'
#' @param features A N x M matrix of item features (Note, this is indeed
#'     a `matrix`, not of class `dist`. This is ensured when calling this
#'     function from the anticlustering wrapper function.)
#' @param distances A N x N matrix representing the
#'     pairwise dissimilarities between all N elements. CAN an be an
#'     object of class \code{dist}.
#' @param K How many anticlusters should be created.
#' @param solver A string identifying the solver to be used ("Rglpk",
#'     "gurobi", or "Rcplex")
#' @param standardize Boolean - should the feature values be
#'     standardized before groups are created?
#' @param preclustering Boolean, should a preclustering be conducted
#'     before anticlusters are created.
#'
#' @return A vector representing the anticluster affiliation of
#'     elements.
#'
#' @noRd

exact_anticlustering <- function(distances, K, solver, preclustering) {

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
