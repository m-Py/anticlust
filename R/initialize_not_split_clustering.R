
# only works for equal-sized groups! (or does it?)
clustering_fit <- function(not_split, cl) {
  not_split <- c(not_split)
  K <- length(unique(cl))
  # each not_split group should be in the same cluster:
  tab <- table(cl, not_split)
  # i.e., in tab, each column should have only zeros except for one row
  (colSums(table(cl, not_split) == 0) - (K - 1))^2
}

obj_fun <- function(x, cl) {
  clustering_fit <- clustering_fit(x, cl)
  if (all(clustering_fit == 0)) {
    return(Inf)
  }
  sum(clustering_fit) * (-1) # clustering fit has to be minimized
}

initialize_with_not_split <- function(not_split, N, K) {
  df <- data.frame(order = 1:N, not_split = not_split)
  df <- df[order(df$not_split), ]
  df$cl <- sort(rep_len(1:K, N))
  if (all(clustering_fit(df$not_split, df$cl) == 0)) {
    print("success")
    return(df[order(df$order), ]$cl)
  }
  # if initial assignment does not work, try to optimize!
  df$cl <- anticlustering(
    df$not_split, 
    K = df$cl, 
    objective = obj_fun,
  )
  if (all(clustering_fit(df$not_split, df$cl) == 0)) {
    print("success after optimization")
  } else {
    print("no success")
  }
  df[order(df$order), ]$cl
}
