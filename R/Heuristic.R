
#' Anticlustering based on a heuristic
#'
#' @param features A data.frame, matrix or vector representing the
#'     features that are used.
#' @param clustering A vector representing the preclustering of
#'     elements.  See Details.
#'
#' @return A vector representing the anticlustering.
#'
#' @details The heuristic approach to anticlustering forbids elements
#'     that are part of the same precluster to be assigned to the same
#'     group. The preclustering should be accomplished by one of the
#'     clustering functions, `equal_sized_cluster_editing` (an exact
#'     method that minimizes distance criterion under the restriction of
#'     equal group sizes) or `equal_sized_clustering` (a heuristic
#'     method that tries to minimize the variance criterion under the
#'     restriction of equal group sizes).
#'
#' @export
#'
heuristic_anticlustering <- function(features, clustering, nrep = 100) {
  ## sort by group, later sort back by item and return group
  dat <- data.frame(group = clustering, features, item = 1:nrow(features))
  dat <- dat[order(dat$group), ]
  ## start optimizing
  best_obj <- -Inf
  n_items <- nrow(dat)
  n_groups <- n_items / length(unique(dat$group))
  ## sequentially try out random assignments that place pregrouped
  ## items in different groups
  for (i in 1:nrep) {
    anticlusters <- unlist(replicate(n_items / n_groups, sample(1:n_groups), simplify = FALSE))
    cur_obj <- obj_value_dist(features, anticlusters)
    if (cur_obj > best_obj) {
      best_assignment <- anticlusters
      best_obj <- cur_obj
    }
  }
  dat$group <- best_assignment
  dat <- dat[order(dat$item), ]
  return(dat$group)
}

#' Objective value for the distance criterion
#'
#' @param features A data.frame, matrix or vector representing the
#'   features that are used in the assignment.
#' @param anticlusters A vector representing the anticluster affiliation
#'
#' @return Scalar, the total sum of within-cluster
#'   pointwise distances.
#'
#' @export
#'

obj_value_distance <- function(features, anticlusters) {
  ## determine distances within each group
  distances <- by(features, assignment, dist)
  ## determine objective as the sum of all distances per group
  objective <- sum(sapply(distances, sum))
  return(objective)
}


#' Heuristic clustering algorithm to create equal-sized clusters
#'
#' Uses kmeans algorithm to initiate cluster centers, then sequentially
#' chooses a cluster center and assigns the closest element to it,
#' ensuring that each cluster is filled with the same number of items.
#'
#' @param features A vector, matrix or data.frame of data points.  Rows
#'     correspond to items and columns correspond to features.
#' @param nclusters The number of clusters to be created. Must be a
#'     factor of the number of items.
#'
#' @return A vector representing the clustering
#'
#' @export
#'
#' @importFrom stats kmeans
#'
equal_sized_kmeans <- function(features, nclusters) {
  if (nclusters <= 1 | nrow(features) %% nclusters != 0)
    stop("The number of features must be a multiplier of nclusters")
  features <- as.matrix(features) #  if features is a vector
  ## initialize cluster centers using kmeans
  centers <- stats::kmeans(features, nclusters)$centers
  ## determine distances between all items and all cluster centers:
  clust_dist <- dist_from_centers(features, centers, squared = TRUE)
  assignments <- heuristic_cluster_assignment(clust_dist)
  ## return in long format that is consistent with return of ILP
  return(clusters_to_long(assignments))
}

## Called from within `equal_sized_kmeans`. Assigns elements to
## clusters, sequentially fills elements to the closest cluster center,
## fills the same number of elements into each cluster
heuristic_cluster_assignment <- function(clust_dist) {
  ## Variables that encode relevant dimensions
  nclusters <- ncol(clust_dist)
  N <- nrow(clust_dist)
  elements_per_cluster <- N / nclusters
  ## Include element ID as column of input because rows are later
  ## removed from matrix
  clust_dist <- cbind(clust_dist, 1:N)
  colnames(clust_dist) <- c(paste0("cluster", 1:nclusters), "item")
  ## Storage of cluster assignments (Columns = clusters; cells indicate the
  ## index of the element that is assigned to the respective cluster)
  assignments <- matrix(nrow = elements_per_cluster, ncol = nclusters)
  for (i in 1:elements_per_cluster) {
    ## iterate over clusters in random order
    for (j in sample(nclusters)) {
      ## 1. Choose element that is closest to center of cluster j:
      min_index <- which.min(clust_dist[, j])
      element <- clust_dist[min_index, "item"]
      ## 2. Add this element to cluster j:
      assignments[i, j] <- element
      ## 3. Remove this element to prevent further assignments:
      clust_dist <- clust_dist[-min_index, , drop = FALSE] # drop = FALSE !
    }
  }
  return(assignments)
}

## Called from within `equal_sized_kmeans`. Turns a matrix of cluster
## assignments into long format.
clusters_to_long <- function(assignments) {
  P <- ncol(assignments) # number of clusters
  N <- P * nrow(assignments)
  N_PER_CLUSTER <- N / P
  long_clusters <- cbind(c(assignments), rep(1:P, each = N_PER_CLUSTER))
  ## Make in correct order:
  long_clusters <- long_clusters[order(long_clusters[, 1]), ]
  return(long_clusters[, 2])
}


#' Edit distances of neighbours
#'
#' Based on a preclustering, distances between items of the same cluster
#' are set to a large negative value. By considering a preclustering of
#' items, very similar items will not be assigned to the same group when
#' the fixed distance object is used to create the ILP formulation of
#' the item assignment instance.
#'
#' @param distances A distance object or matrix of
#'     between-item-distances.
#' @return A vector representing the group assignment; objects that are
#'     part of the same group as indicated by this vector are assigned a
#'     new distance.
#' @param value The value that is assigned to the fixed distances.
#'     Defaults to -1,000,000 currently.
#'
#' @return A distance object containing the fixed distances. For items
#'     that are part of the same cluster (as specified in the
#'     `data.frame` `assignment`), the between-item-distances are set to
#'     -1000000. This will have to be replaced by a theoretically-sound
#'     value; -1000000 is just a hack that will work in the present
#'     applications.
#'
#' @importFrom stats as.dist dist
#' @export

edit_distances <- function(distances, assignment, value = -1000000) {
  n_groups <- length(unique(assignment))
  distances <- as.matrix(distances)
  for (i in 1:n_groups) {
    items <- which(assignment == i) ## which items are in the i'th cluster
    ## two for-loops to tap all pairwise distances that may require
    ## fixing (a lot of unnecessary iterations, probably)
    for (j in 1:(length(items) - 1)) {
      for (t in 2:length(items)) {
        distances[items[j], items[t]] <- value
        distances[items[t], items[j]] <- value
      }
    }
  }
  return(as.dist(distances))
}
