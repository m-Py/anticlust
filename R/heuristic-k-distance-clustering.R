
greedy_k_clustering <- function(distances, K) {
  distances <- as.matrix(distances)
  diag(distances) <- Inf
  N <- nrow(distances)

  ## Columns = anticlusters, rows = elements
  clusters <- init_clusters(distances, N, K)
  ## remove used stimuli from the pool
  complete <- clusters[complete.cases(clusters), ]
  up_to_now <- as.matrix(expand.grid(c(complete), c(complete)))
  distances[up_to_now] <- Inf

  ## fill the rest
  for (i in 3:(N/K)) {
    for (k in 1:K) {
      ## Select items that are already in k'th cluster
      items_in_cluster <- clusters[, k]
      items_in_cluster <- items_in_cluster[!is.na(items_in_cluster)]
      ## Which is the new item that is inserted:
      insert <- find_closest_to_cluster(distances, items_in_cluster)
      ## append new item to k'th anticluster
      clusters[i, k] <- insert
      ## Remove item from pool:
      complete <- c(clusters[!is.na(clusters)])
      up_to_now <- as.matrix(expand.grid(insert, c(complete)))
      distances[up_to_now] <- Inf
      distances[up_to_now[, 2:1]] <- Inf
    }
  }

  ## Convert matrix of clusters into vector
  groups <- rep(NA, N)
  for (k in 1:K) {
    groups[clusters[, k]] <- k
  }

  ## Assert that the output has the correct structure:
  if (any(table(groups) != table(groups)[2])) {
    stop("Something went wrong: not the same anticluster sizes")
  }
  if (any(is.na(groups))) {
    stop("Something went wrong: NA in output")
  }

  return(groups)
}


#' For a cluster given items that are already in the cluster,
#' which other item should be included. Currently this is *very*
#' greedy, only looking for the minimum distance between two elements.
#' But I want the minimum distance between all elements of the cluster
#' and another element
find_closest_to_cluster <- function(distances, items_in_cluster) {
  ## select items that are connected to the items in the current cluster
  current_min <- min(distances[, items_in_cluster])
  ## select only relevant columns
  tmp_distances <- distances
  tmp_distances[, -items_in_cluster] <- Inf
  insert <- which(tmp_distances == current_min, arr.ind = TRUE)
  insert <- setdiff(insert, items_in_cluster)
  return(insert)
}



## Initialize assignment: Fill K anticlusters with two elements
init_clusters <- function(distances, N, K) {
  anticlusters <- matrix(nrow = N / K, ncol = K)
  ## Pick the first two elements per anticluster
  for (k in 1:K) {
    # select the first pair if a duplicate of the minimum occurs:
    pair <- which(distances == min(distances), arr.ind = TRUE)[1, ]
    anticlusters[1:2, k] <- pair
    distances[pair, ] <- Inf
    distances[, pair] <- Inf
  }
  return(anticlusters)
}

# Create artifical data
N <- 21
K <- 3
n_features <- 2
features <- matrix(rnorm(N * n_features), ncol = n_features)

tt <- as.matrix(dist(features))

library("anticlust")

ac_greedy <- greedy_k_clustering(tt, K)
ac_var <- balanced_clustering(features, K = K, method = "heuristic")
obj_shit <- anticlust:::obj_value_distance(features, ac_greedy)
obj_var <- anticlust:::obj_value_distance(features, ac_var)
ac_exact <- balanced_clustering(features, K = K, method = "exact")
obj_good <- anticlust:::obj_value_distance(features, ac_exact)

#ac_greedy <- (ac_greedy - 1); ac_greedy[ac_greedy == 0] <- 2
par(mfrow = c(1, 2))
plot_clusters(features, ac_greedy, cex = 2, main = "greedy", xlab = "", ylab = "")
plot_clusters(features, ac_exact, cex = 2, main = "exact", xlab = "", ylab = "")

obj_good / obj_shit
obj_good / obj_var


#' Determine the minimum (or maximum) sum of distances between
#' elements in a cluster and elements that are not yet in a cluster
#'
#' @param distances A distance matrix. It is in preprocessed form so
#'     that several entries will be Inf or -Inf.
#' @param cluster_elements The elements that are in the cluster under
#'     processing.
#'
#' @return The element that is closest to (or furthest apart from)
#'     a cluster.
#'
