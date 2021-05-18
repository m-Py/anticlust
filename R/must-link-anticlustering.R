
# Quantify how well a clustering satisfies must-link constraints
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
  # clustering fit has to be minimized, but my exchange method maximizes, 
  # so return the objective function (*-1)
  sum(clustering_fit) * (-1)
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

# Given an assignment of original elements to clusters, generate the
# corresponding clusters for the merged elements. 
original_cluster_to_merged_cluster <- function(clusters, must_link) {
  df <- data.frame(clusters, must_link)
  new_df <- df[!duplicated(df$must_link), ]
  # guarantee order by must_link grouping! 
  # (probably after merge it already is sorted this way, though)
  new_df[order(new_df$must_link), ]$clusters
}

# Given an assignment of "merged" elements to clusters, re-establish the
# corresponding clusters for the original elements. 
merged_cluster_to_original_cluster <- function(merged_clusters, must_link) {
  df <- data.frame(must_link, order = 1:length(must_link))
  new_order <- order(must_link)
  df <- df[new_order, ]
  df$clusters <- rep(merged_clusters, table(must_link))
  df[order(df$order), "clusters"]
}

#' @export
anticlustering_notsplit <- function(x, K, must_link) {
  
  clusters_init <- initialize_must_link_clustering(must_link, N = NROW(x), K = K)
  
  DF_ <- data.frame(must_link, x)
  
  obj_for_merged_clusters <- function(x, clusters) {
    clusters_real <- merged_cluster_to_original_cluster(clusters, DF_[, 1])
    # only accept legal exchanges!
    if (any(sort(table(clusters_real)) != sort(table(clusters_init)))) {
      return(-Inf)
    }
    kplus_objective(DF_[, -1], clusters_real)
  }
  
  dummy_data <- 1:length(unique(must_link))
  dummy_groups <- anticlustering(
    dummy_data,
    K = original_cluster_to_merged_cluster(clusters_init, must_link),
    objective = obj_for_merged_clusters
  )
  merged_cluster_to_original_cluster(dummy_groups, must_link)
}

