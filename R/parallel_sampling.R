
parallel_sampling <- function(dat, K, objective, nrep, sampling_plan,
                              use_distances) {
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  reps_per_cluster <- ceiling(nrep / ncores)
  ## this is missing apparently:
  # parallel::clusterExport(
  #   cl = cl,
  #   varlist = c("..."),
  #   envir = environment()
  # )


  assignments <- parallel::parLapply(
    X = 1:ncores,
    fun = random_sampling,
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
    fun = obj_value,
    data = dat[, -(1:2), drop = FALSE]
  )
  best_obj <- which.max(objectives)
  assignments[[best_obj]]
}
