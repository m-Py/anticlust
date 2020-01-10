
# Internal function for subset selection based on preclustering - then 
# anticlustering. if `split_by` exists, min-max anticlustering is conducted
subset_anticlustering <- function(data, split_by, equalize, balance, design, n) {
  
  # total N of input data:
  N <- nrow(data)
  # number of groups that are needed:
  K <- prod(design)
  
  # track which items are selected via ID:
  data$id_anticlustering <- 1:N
  
  # test if I can safely remove all entries with NA values
  preselection <- safely_exclude_na(data, split_by, equalize, balance, K, n)
  
  # get most similar items, possibly under categorical restrictions
  preclusters <- categorical_restrictions(preselection, equalize, balance, K)
  
  # only select clusters of size `K`
  filled <- which(table(preclusters) == K)
  is_selected <- preclusters %in% filled
  preselection <- preselection[is_selected, , drop = FALSE]
  preclusters <- preclusters[is_selected]
  
  # Get best clusters - similar wrt `equalize`; dissimilar wrt `split_by`
  distances <- by(preselection[, equalize], preclusters, dist)
  objectives <- sapply(distances, sum)
  
  # Some additional work needs to be done if min-max anticlustering is required
  if (argument_exists(split_by)) {
    distances2 <- by(preselection[, split_by], preclusters, dist)
    objectives2 <- sapply(distances2, sum)
    objectives <- objectives - objectives2
  }
  
  needed_n <- K * n
  needed_clusters <- needed_n / table(preclusters)[1]
  # Best preclusters
  cluster_ids <- 1:needed_clusters
  most_similar_clusters <- as.numeric(names(sort(objectives))[cluster_ids])
  
  # encode items that are selected (all that are in the best preclusters)
  is_in_output <- preclusters %in% most_similar_clusters
  preselection <- preselection[is_in_output, ]
  preclusters  <- preclusters[is_in_output]
  
  iv <- NULL
  obj_fun <- "distance"
  if (argument_exists(split_by)) {
    iv <- scale(preselection[, split_by])
    obj_fun <- make_obj_function(data, equalize, split_by, design) 
  }
  
  # now optimize assignment using exchange method
  anticlusters <- anticlustering(
    features = scale(preselection[, c(equalize, split_by)]),
    K = K,
    categories = preclusters,
    objective = obj_fun
  )
  
  # prepare output: Needs to incorporate NA for non-selected items
  output <- rep(NA, N)
  output[preselection$id_anticlustering] <- anticlusters
  output
}

safely_exclude_na <- function(data, split_by, equalize, balance, K, n) {
  has_no_na <- complete.cases(data[, equalize, drop = FALSE])
  if (argument_exists(split_by)) {
    has_no_na <- has_no_na & complete.cases(data[, split_by, drop = FALSE])
  }
  if (argument_exists(balance)) {
    has_no_na <- has_no_na & complete.cases(data[, balance, drop = FALSE])
  }
  # can remove data if enough entries remain after NA exclusion
  if ((K * n) <= sum(has_no_na)) {
    data <- data[has_no_na, , drop = FALSE]
    if (sum(!has_no_na) > 0) {
      message("\n", sum(!has_no_na), " records could be excluded due to missing values.\n",
              "Selecting stimuli from the remaining ", nrow(data), " records.")
    }
  }
  data
}


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
