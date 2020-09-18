
# Function that repeatedly calls anticlustering and returns best results
repeat_anticlustering <- function(x, K, objective, categories, preclustering, 
                                  method, repetitions) {
  
  if (inherits(objective, "function")) {
    obj_function <- objective
  } else if (objective == "variance") {
    obj_function <- variance_objective
  } else if (objective == "diversity" || objective == "distance") {
    obj_function <- diversity_objective
  } else if (objective == "dispersion") {
    obj_function <- dispersion_objective
  }
  
  if (method == "local-maximum") {
    candidate_solutions <- lapply(
      1:repetitions, 
      local_maximum_anticlustering,
      data = x, 
      K = K, 
      objective = objective, 
      obj_function = obj_function,
      categories = categories, 
      preclustering = preclustering
    )
  } else if (method == "exchange") {
    candidate_solutions <- lapply(
      1:repetitions, 
      anticlustering_,
      data = x, 
      K = K, 
      objective = objective, 
      categories = categories, 
      preclustering = preclustering
    )
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


# Function that calls anticlustering exchange method until a local 
# maximum is reached (i.e., no possible swap improves the objective)
local_maximum_anticlustering <- function(X, data, K, objective, obj_function, 
                                         categories, preclustering) {
  
  exchange_partners <- get_categorical_constraints(data, K, preclustering, categories)
  clusters <- initialize_clusters(nrow(data), K, exchange_partners)

  new_obj <- obj_function(data, clusters)
  old_obj <- -Inf
  
  while (new_obj > old_obj) {
    clusters <- anticlustering(
      x = data, 
      K = clusters, 
      objective = objective, 
      categories = exchange_partners, 
      preclustering = FALSE,
      repetitions = NULL,
      method = "exchange"
    )
    old_obj <- new_obj
    new_obj <- obj_function(data, clusters)
  }
  clusters
}

# `anticlustering()` to be used repeatedly by `lapply()`
anticlustering_ <- function(X, data, K, objective, preclustering, categories) {
  anticlustering(
    x = data, 
    K = K, 
    objective = objective, 
    method = "exchange", 
    preclustering = preclustering, 
    categories = categories, 
    repetitions = NULL
  )
}
