
library("anticlust")

# exported function cannot be used with NA"
data(iris)
iris[1, 1] <- NA
expect_error(
  anticlustering(iris[, -5], K = 3)
)
expect_error(
  matching(iris[, -5], p = 3)
)
expect_error(
  balanced_clustering(iris[, -5], p = 3)
)

# exported functions cannot be used with non-numeric input"
data(iris)
expect_error(
  anticlustering(iris, K = 3)
)
expect_error(
  matching(iris, p = 3)
)
expect_error(
  balanced_clustering(iris, p = 3)
)
