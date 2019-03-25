
# Return group membership based on ILP solution for item assignment
#
# This function is primarily called from other functions and
# probably should not be called directly. Only use this if you have
# called solve_ilp earlier.
#
# @param ilp An ILP representation of the item assignment problem
#     (returned by anticlustering_ilp)
# @param solution The solution of the instance returned by
#     solve_ilp
#
# @return A vector representing anticluster affiliation for each
#     element
#
ilp_to_groups <- function(ilp, solution) {
  assignment_list <- fill_groups_wce_(ilp, solution$x)
  ret <- data.frame(item = unlist(assignment_list),
                    group = rep(1:length(assignment_list),
                                list_length(assignment_list)))
  ret <- sort_by_col(ret, "item")
  return(ret$group)
}


## Extract groups from a the information of pairwise connectivity of
## items
fill_groups_wce_ <- function(ilp, variable_assignment) {
  ## as long as this function is not in the package: I have to load
  ## the functions in file ILP.R manually first
  n_items <- length(union(ilp$costs$i, ilp$costs$j))
  groups <- as.list(1:n_items) # initialize groups: each item is a "cluster leader"
  groups <- fill_groups_(ilp, variable_assignment, groups)

  ## Clean output of `fill_groups_`, which assumes (a) cluster leaders
  ## and (b) equally sized clusters.

  ## First step: Assign items to empty elements

  for (i in 1:length(groups)) {
    if (length(groups[[i]]) == 0)
      groups[[i]] <- i ## if i: item is its own cluster; if 0: remove cluster
  }

  ## Second step: Mark all duplicated clusters

  for (i in 1:length(groups)) {
    ## this has to be redone each iteration to ensure that the last
    ## duplicate is not marked (it has to remain!):
    duplicates <- duplicated(groups)
    current_marker <- duplicates[i]
    if (current_marker == TRUE) { # is duplicate
      groups[[i]] <- 0 # mark with 0
    }
  }

  ## Third step: Remove marked duplicates
  marked <- sapply(groups, function(x) identical(x, 0))
  groups <- groups[!marked]

  return(groups)
}

## This finds which items are in the same cluster, but it possibly
## returns duplicates
fill_groups_ <- function(ilp, variable_assignment, groups) {
  costs <- ilp$costs
  for (i in 1:length(groups)) {
    leader <- groups[[i]][1]
    ## select all edges of the "leader"
    edges  <- costs[(costs$i == leader | costs$j == leader) & variable_assignment == 1, ]
    group_elements <- unique(c(edges$i, edges$j))
    groups[[i]] <- group_elements
  }
  ## return as data.frame
  return(groups)
}

## get the length of each vector element of a list:
list_length <- function(x) sapply(x, length)


# Edit distances
#
# Based on a clustering, distances between elements of the same cluster
# are set to a specified value. This is used to enforce or forbid pairs
# of elements to be part of the same anticluster when the distance
# matrix is passed to the ILP solver.  By considering a preclustering
# of items, very similar items will not be assigned to the same group
# when the fixed distance object is used to create the ILP formulation
# of the item assignment instance.
#
# @param distances A n xn matrix of between-item-distances.
#
# @param clustering A vector representing a clustering; objects that
#     are part of the same cluster are assigned a new distance.
# @param value The value that is assigned to the fixed distances.
#     Defaults to -1,000,000.
#
# @return A distance object containing the fixed distances. For items
#     that are part of the same cluster.
#
#
edit_distances <- function(distances, clustering, value = -1000000) {
  n_groups <- length(unique(clustering))
  for (i in 1:n_groups) {
    items <- which(clustering == i) ## which items are in the i'th cluster
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
