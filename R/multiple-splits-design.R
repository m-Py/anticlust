

# Generate an objective functions for min-max anticlustering

make_obj_function <- function(data, equalize, split_by, design) {
  
  # get all levels of clusters and the corresponding levels for split variables
  levels <- expand.grid(lapply(design, function(x) 1:x))
  
  # combine to a single objective function
  function(cl, data) {
    similarity_covariates <- obj_value_distance(
      cl, 
      data[, equalize, drop = FALSE]
    )
    dissims <- c()
    for (i in 1:length(design)) {
      clusters <- levels[, i][cl]
      # compute sum of distances within cluster
      distances <- by(data[, split_by[i]], clusters, dist)
      dissims[i] <- sum(sapply(distances, sum))
    }
    similarity_covariates - sum(dissims)
  }
}
