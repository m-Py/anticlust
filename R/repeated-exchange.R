
# Function that repeatedly calls anticlustering and returns best results
repeat_anticlustering <- function(x, K, objective, categories, method, repetitions) {
  
  N <- nrow(x)
  
  # Create initial cluster assignment for each `repetition`
  if (argument_exists(repetitions) && repetitions > 1) {
    clusters <- get_multiple_initial_clusters(N, K, categories, repetitions)
  } else {
    clusters <- list(initialize_clusters(N, K, categories))
  }
  
  if (inherits(objective, "function")) {
    obj_function <- objective
  } else if (objective == "variance") {
    obj_function <- variance_objective
  } else if (objective == "diversity" || objective == "distance") {
    obj_function <- diversity_objective
  } else if (objective == "dispersion") {
    obj_function <- dispersion_objective
  } else if (objective == "average-diversity") {
    obj_function <- function(clusters, x) {
      weighted_diversity_objective_(clusters, x, table(clusters))
    }
  }
  
  if (method == "local-maximum") {
    candidate_solutions <- lapply(
      clusters,
      local_maximum_anticlustering,
      data = x, 
      objective = objective, 
      obj_function = obj_function,
      categories = categories
    )
  } else if (method == "exchange") {
    candidate_solutions <- lapply(
      clusters,
      anticlustering_,
      data = x, 
      objective = objective, 
      categories = categories
    )
    candidate_solutions
  }
  
  # Define objective function that can be used by `lapply` on many clusterings
  obj_function_ <- function(clusters, x) {
    obj_function(x, clusters)
  }
  # Get best of all solutions
  objs <- lapply(
    candidate_solutions,
    obj_function_,
    x = x
  )
  candidate_solutions[[which.max(objs)]]
}

get_multiple_initial_clusters <- function(N, K, categories, repetitions) {
  
  clusters <- list(initialize_clusters(N, K, categories))
  
  if (argument_exists(categories)) {
    clusters_repetitions <- as.list(data.frame(as.matrix(
      replicate(repetitions - 1, categorical_sampling(categories, K))
    )))
  } else {
    clusters_repetitions <- as.list(data.frame(as.matrix(
      replicate(repetitions - 1, sample(clusters[[1]]))
    )))
  }
    
  clusters_repetitions <- unname(clusters_repetitions)
  merge_lists(list(clusters, clusters_repetitions))
}


# Function that calls anticlustering exchange method until a local 
# maximum is reached (i.e., no possible swap improves the objective)
local_maximum_anticlustering <- function(
  clusters, data, objective, obj_function, categories) {
  
  old_obj <- obj_function(data, clusters)
  has_improved <- TRUE
  
  while (has_improved) {
    new_clusters <- anticlustering(
      x = data, 
      K = clusters, 
      objective = objective, 
      categories = categories, 
      repetitions = NULL,
      method = "exchange",
      preclustering = FALSE
    )
    new_obj <- obj_function(data, new_clusters)
    if (new_obj > old_obj) {
      clusters <- new_clusters
      old_obj <- new_obj
    } else {
      has_improved <- FALSE
    }
  }
  clusters
}

# `anticlustering()` to be used repeatedly by `lapply()`
anticlustering_ <- function(clusters, data, objective, categories) {
  anticlustering(
    data, 
    K = clusters, 
    objective = objective, 
    method = "exchange", 
    categories = categories, 
    repetitions = NULL,
    preclustering = FALSE
  )
}
