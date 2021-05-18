
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
