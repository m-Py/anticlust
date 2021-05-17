
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

