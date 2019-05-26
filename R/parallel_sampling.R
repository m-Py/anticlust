
parallel_sampling <- function(dat, K, objective, nrep, sampling_plan,
                              use_distances) {
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  reps_per_cluster <- ceiling(nrep / ncores)

  assignments <- parallel::parLapply(
    X = 1:ncores,
    fun = lapply_random_samling,
    cl = cl,
    dat = dat,
    K = K,
    objective = dat,
    nrep = reps_per_cluster,
    sampling_plan = sampling_plan,
    use_distances = use_distances
  )
  parallel::stopCluster(cl)

  if (objective == "distance" && use_distances == TRUE) {
    obj_value <- distance_objective_
  } else if (objective == "distance" && use_distances == FALSE) {
    obj_value <- obj_value_distance
  } else {
    obj_value <- variance_objective_
  }

  objectives <- lapply(
    assignments,
    FUN = obj_value,
    data = dat[, -(1:2), drop = FALSE]
  )
  best_obj <- which.max(objectives)
  assignments[[best_obj]]
}

# Make random_sampling usable for lapply (additional argument x)
lapply_random_samling <- function(x, dat, K, objective, nrep, sampling_plan,
                                use_distances) {
  random_sampling(dat, K, objective, nrep, sampling_plan,
                  use_distances)
}

