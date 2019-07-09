
#' Parallel random sampling for anticlustering
#'
#' A wrapper function that sets up some clusters and sends part of
#' all required reptitions to each cluster. In the end, the best solution
#' across the solutions found by each cluster is returned.
#'
#' @inheritParams random_sampling
#' @param seed A seed passed to \code{parallel::clusterSetRNGStream}
#'
#' @return The best anticlustering solution
#'
#' @noRd
#'
#' @importFrom parallel detectCores makeCluster clusterSetRNGStream parLapply stopCluster
#'

parallel_sampling <- function(dat, K, nrep, sampling_plan,
                              obj_function, seed, ncores = NULL) {
  if (!argument_exists(ncores)) {
    ncores <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(ncores)
  reps_per_cluster <- ceiling(nrep / ncores)

  if (argument_exists(seed)) {
    clusterSetRNGStream(cl, seed)
  }

  assignments <- parallel::parLapply(
    X = 1:ncores,
    fun = lapply_random_samling,
    cl = cl,
    dat = dat,
    K = K,
    nrep = reps_per_cluster,
    sampling_plan = sampling_plan,
    obj_function = obj_function
  )
  on.exit(parallel::stopCluster(cl))

  ## Compute objectives for the best values of each core:
  objectives <- lapply(
    assignments,
    FUN = obj_function,
    data = dat[, -(1:2), drop = FALSE]
  )
  best_obj <- which.max(objectives)
  assignments[[best_obj]]
}

# Make random_sampling usable for lapply (additional argument x)
lapply_random_samling <- function(x, dat, K, nrep, sampling_plan,
                                  obj_function) {
  random_sampling_(dat, K, nrep, sampling_plan,
                  obj_function)
}

