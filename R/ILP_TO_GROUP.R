
#######

## The following functions take a solution of the item assignment
## problem and return an assignment of items to groups as a data.frame.
## (function `ilp_to_groups` is the highest level function that should
## be called)

#' Return group membership based on ILP solution for an item assignment
#' instance
#'
#' This function is primarily called from within the function
#' `item_assignment`.
#'
#' @param ilp An ILP representation of the item assignment problem
#'     (returned by rom `item_assign_ilp`)
#' @param solution The solution of the instance returned by the function
#'     `solve_ilp`.
#' @param items The data.frame that was originally passed to
#'     item_assignment containing all item features.
#'
#' @return A data.frame containing one column of item ids (here, the id
#'     corresponds to the order of the items) and one column contains the
#'     group assignments of the items. The original items are also
#'     returned as columns of the data.frame.
#' @export
#'

ilp_to_groups <- function(ilp, solution, items) {
  assignment <- fill_groups_wce(ilp, solution)
  ## the following works because correct order of items is ensured
  ## through function group_assignment, that sorts by item
  assignment <- data.frame(assignment, items)
  return(assignment)
}


## TODO: better name and comment the functions below; they replace
## several previous functions that were made obsolete by the more
## general functions

#' Return group membership based on ILP solution for a weighted cluster
#' editing instance
#'
#' @param ilp An ILP representation of the item assignment problem
#'     (returned by rom `item_assign_ilp`)
#' @param solution The solution of the instance returned by the function
#'     `solve_ilp`.
#'
#' @return A data.frame containing one column of item ids (here, the id
#'     corresponds to the order of the items); one column contains the
#'     group assignments of the items.
#' @export
#'
fill_groups_wce <- function(ilp, solution) {
  assignment_list <- fill_groups_wce_(ilp, solution$x)
  ret <- data.frame(item = unlist(assignment_list),
                    group = rep(1:length(assignment_list),
                                list_length(assignment_list)))
  ret <- ret[order(ret$item), ]
  rownames(ret) <- NULL
  return(ret)
}


## extract groups from a weighted cluster editing instance (i.e., there
## are no cluster leaders an groups may differ in length). It actually calls
## `fill_groups_` and cleans the output.
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

## This finds which items are in the same cluster, but
## it possibly returns duplicates
fill_groups_ <- function(ilp, variable_assignment, groups) {
  select <- grepl("x", names(variable_assignment)) ## disregard cluster leader variables
  ## connections = assignment of all x_ij decision variables
  connections <- variable_assignment[select]
  ## check out what items each cluster leader is connected to
  for (i in 1:length(groups)) {
    leader <- groups[[i]][1]
    ## select all edges of the leader:
    edges  <- subset(ilp$costs, (i == leader | j == leader) & connections == 1)
    group_elements <- unique(c(edges$i, edges$j))
    groups[[i]] <- group_elements
  }
  ## return as data.frame
  return(groups)
}

## get the length of each vector element of a list:
list_length <- function(x) sapply(x, length)
