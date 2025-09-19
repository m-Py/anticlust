
library("anticlust")

# some exported function cannot be used with NA (anticlustering now can!)
data(iris)
iris[1, 1] <- NA
expect_error(
  matching(iris[, -5], p = 3)
)
expect_error(
  balanced_clustering(iris[, -5], K = 3)
)

# some exported functions cannot be used with non-numeric input (anticlustering now can!)
data(iris)
expect_error(
  matching(iris, p = 3)
)
expect_error(
  balanced_clustering(iris, K = 3)
)
