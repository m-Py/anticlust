
#' @export
cannot_link_anticlustering <- function(x, init_clusters, cannot_link, objective, ...) {

  if (objective == "kplus") {
    x <- kplus_moment_variables(x, 2)
    objective <- "variance"
  }
  
  if (objective == "variance") {
    x <- convert_to_distances(x)^2
  } else if (objective == "diversity") {
    x <- convert_to_distances(x)
  } else {
    stop("Argument 'objective' must be 'diversity', 'variance' or 'kplus'")
  }
  frequencies <- table(init_clusters)
  if (!all(frequencies == frequencies[1])) {
    stop("Currently, cannot-link anticlustering is only supported for equal-sized groups.")
  }
  

  # set cannot-link distances to large negative value so they cannot be linked
  x[rbind(cannot_link, t(apply(cannot_link, 1, rev)))] <- -(sum(x) + 1)
  
  
  ## special case of only one init partition...
  init_clusters <- as.matrix(init_clusters)
  if (ncol(init_clusters) == 1) {
    init_clusters <- t(init_clusters)
  }
  
  solutions <- list()
  objectives <- rep(NA, nrow(init_clusters))
  
  for (i in 1:nrow(init_clusters)) {
    solutions[[i]] <- anticlustering(x, K = init_clusters[i, ], objective = "diversity", ...)
    objectives[i] <- diversity_objective_(solutions[[i]], x)
  }
  solutions[[which.max(objectives)]]
  
}
