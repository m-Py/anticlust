
# only works for equal-sized groups! (or does it?)
clustering_fit_must_link <- function(must_link, cl) {
  must_link <- c(must_link)
  K <- length(unique(cl))
  # each not_split group should be in the same cluster
  # (# i.e., in tab, each column should have only zeros except for one row):
  tab <- table(cl, must_link)
  # compute squared difference between the number of zeros in each column
  # and the number of zeros that should be there; this is the objective function
  # that has to be minimized
  (colSums(table(cl, must_link) == 0) - (K - 1))^2
}

# Expand the above objective function, return Infinity when the clustering
# fits perfectly.
obj_fun_must_link <- function(x, cl) {
  clustering_fit <- clustering_fit_must_link(x, cl)
  if (all(clustering_fit == 0)) {
    return(Inf)
  }
  sum(clustering_fit) * (-1) # clustering fit has to be minimized
}

# Initialize a clustering vector satisfying must-link constraints.
# If an initial guess does not work, uses an exchange process to optimize 
# the above objective function representing the degree to which the constraints
# are satisfied.
initialize_must_link_clustering <- function(must_link, N, K) {
  df <- data.frame(order = 1:N, must_link = must_link)
  df <- df[order(df$must_link), ]
  df$cl <- sort(rep_len(1:K, N))
  if (all(clustering_fit_must_link(df$must_link, df$cl) == 0)) {
    print("success")
    return(df[order(df$order), ]$cl)
  }
  # if initial assignment does not work, try to optimize!
  df$cl <- exchange_method(
    as.matrix(df$must_link), 
    K = df$cl, 
    obj_function = obj_fun_must_link,
    categories = NULL
  )
  if (all(clustering_fit_must_link(df$must_link, df$cl) == 0)) {
    print("success after optimization")
  } else {
    print("no success")
  }
  df[order(df$order), ]$cl
}


# argument `clusters` is an assignment based on matches
merged_cluster_to_original_cluster <- function(merged_clusters, not_split) {
  df <- data.frame(not_split, order = 1:length(not_split))
  new_order <- order(not_split)
  df <- df[new_order, ]
  df$clusters <- rep(merged_clusters, table(not_split))
  df[order(df$order), "clusters"]
}

#' @export
anticlustering_notsplit <- function(x, K, not_split) {
  
  clusters_init <- initialize_clusters(NROW(x), K = K, NULL)
  
  DF_ <- data.frame(not_split, x)
  
  obj_for_merged_clusters <- function(x, clusters) {
    clusters_real <- merged_cluster_to_original_cluster(clusters, DF_[, 1])
    # only accept legal exchanges!
    if (any(sort(table(clusters_real)) != sort(table(clusters_init)))) {
      print(table(clusters_real))
      print(table(clusters_init))
      return(-Inf)
    }
    kplus_objective(DF_[, -1], clusters_real)
  }
  
  dummy_data <- 1:length(unique(not_split))
  dummy_groups <- anticlustering(
    dummy_data,
    K = clusters_init,
    objective = obj_for_merged_clusters
  )
  merged_cluster_to_original_cluster(dummy_groups, not_split)
}

