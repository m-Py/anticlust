
# Initialize must-link constraints by assigning all persons having the same ID to the same set
# all others remain free (i.e., as NA)
get_init_assignments <- function(N, ID, target_groups) {
  # Initialize all as NA
  init <- rep(NA, N)
  # for loop to satisfy must-link constraints, will not work in all cases (this can
  # be solved optimally in polynomial time but I do not have this algorithm right now):
  K <- length(target_groups)
  cluster_sizes_real <- rep(0, K)
  multiple_IDs <- as.numeric(names(table(ID)[table(ID) > 1]))
  
  for (current_id in multiple_IDs) {
    random_order_clusters <- sample(K)
    for (k in random_order_clusters) {
      # only fill into cluster if it fits
      if ((cluster_sizes_real[k] + sum(ID == current_id)) > target_groups[k]) {
        if (k == random_order_clusters[K]) {
          stop("Sorry, this failed")
        }
        next
      }
      init[ID == current_id] <- k
      cluster_sizes_real[k] <- cluster_sizes_real[k] + sum(ID == current_id)
      break
    }
  }
  stopifnot(sum(!is.na(init))  == sum(ID %in% multiple_IDs))
  init
}

# After initial assignment, fill the rest randomly
fill_groups <- function(init, target_groups) {
  K <- length(target_groups)
  table_assigned <- table(init)
  # assign elements that have no group (unfortunately, this "simple" task is quite difficult in general)
  if (length(table_assigned) != length(target_groups)) {
    table_assigned <- data.frame(K = 1:K, size = 0)
    df <- as.data.frame(table(init))
    together <- merge(table_assigned, df, by.x = "K", by.y = "init", all = TRUE)
    table_assigned <- ifelse(is.na(together$Freq), 0, together$Freq)
  }
  freq_not_assigned <- target_groups - table_assigned
  init[is.na(init)] <- anticlust:::sample_(rep(1:K, freq_not_assigned))
  stopifnot(all(table(init) == target_groups))
  init
}

# wrapper for the two above
init_must_link_groups <- function(N, ID, target_groups) {
  init <- get_init_assignments(N, ID, target_groups)
  fill_groups(init, target_groups)
}

# Given an assignment of original elements to clusters, generate the
# corresponding clusters for the merged elements. 
original_cluster_to_merged_cluster <- function(clusters, must_link) {
  df <- data.frame(clusters, must_link)
  new_df <- df[!duplicated(df$must_link), ]
  # order by must_link grouping
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

#' Must-link anticlustering
#' 
#' @param x A distance matrix (Only maximum diversity is implemented)
#' @param K What do the required groups look like (as in \code{\link{anticlustering}})
#' @param must_link Vector of IDs; elements having the same value are linked in the same anticluster.
#' @param method  "exchange" (default) or "local-maximum"
#' @export
#' 
#' @details Conducts an exchange method on "merged" elements.
#' 
must_link_anticlustering <- function(x, K, must_link, method = "exchange") {
  validate_input(method, "method", input_set = c("local-maximum", "exchange"), not_na = TRUE, len = 1) 
  stopifnot(is_distance_matrix(x))
  x <- to_matrix(x)
  N <- nrow(x)
  
  clusters_init <- init_must_link_groups(N, must_link, table(initialize_clusters(N, K, NULL)))

  DF_ <- data.frame(must_link, x)
  
  obj_for_merged_clusters <- function(x, clusters) {
    clusters_real <- merged_cluster_to_original_cluster(clusters, DF_[, 1])
    cluster_sizes_feasible <- same_cluster_sizes(clusters_real, clusters_init)
    if (!cluster_sizes_feasible) {
      return(-Inf)
    } 
    diversity_objective(DF_[, -1], clusters_real)
  }
  
  dummy_data <- 1:length(unique(must_link))
  dummy_groups <- anticlustering(
    dummy_data,
    K = original_cluster_to_merged_cluster(clusters_init, must_link),
    objective = obj_for_merged_clusters, 
    method = method
  )
  solution <- merged_cluster_to_original_cluster(dummy_groups, must_link)
  if (!is.infinite(obj_fun_must_link(must_link, solution))) {
    stop("I could not fulfil the `must_link` restrictions, sorry!")
  }
  # Group sizes may be unequal if restrictions were not fulfilled
  if (!same_cluster_sizes(solution, clusters_init)) {
    stop("I could not fulfil the `must_link` restrictions, sorry!")
  }
  solution
}

same_cluster_sizes <- function(clusters1, clusters2) {
  all(sort(table(clusters1)) == sort(table(clusters2)))
}




############################# DEAD CODE:

# Quantify how well a clustering satisfies must-link constraints
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
initialize_must_link_clustering <- function(must_link, N, K, categories) {
  df <- data.frame(order = 1:N, must_link = must_link)
  if (argument_exists(categories)) {
    df$categories <- categories
  }
  df <- df[order(df$must_link), ]
  df$cl <- initialize_clusters(N, K, df$categories)
  if (all(clustering_fit_must_link(df$must_link, df$cl) == 0)) {
    print("success")
    return(df[order(df$order), ]$cl)
  }
  # if initial assignment does not work, try to optimize!
  df$cl <- exchange_method(
    as.matrix(df$must_link), 
    K = df$cl, 
    obj_function = obj_fun_must_link,
    categories = df$categories
  )
  if (all(clustering_fit_must_link(df$must_link, df$cl) == 0)) {
    print("success after optimization")
  } else {
    print("no success")
  }
  df[order(df$order), ]$cl
}


########################## END DEAD