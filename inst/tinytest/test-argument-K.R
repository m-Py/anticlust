
library("anticlust")

# Test standard usage: K = number of groups
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = 3,
  categories = schaper2019$room
)
expect_true(all(table(anticlusters) == c(32, 32, 32)))

expect_true(
  all(table(anticlusters, schaper2019$room) == matrix(rep(16, 6), ncol = 2))
)

# Test K = group size
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = c(48, 24, 24),
  categories = schaper2019$room
)
expect_true(all(table(anticlusters) == c(48, 24, 24)))

# somewhat more extreme values:
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = c(80, 8, 8),
  categories = schaper2019$room
)

expect_true(
  all(table(anticlusters, schaper2019$room) == matrix(c(40, 4, 4, 40, 4, 4), ncol = 2))
)

# Test K = a priori assignment of cases
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = anticlusters,
  categories = schaper2019$room
)
expect_true(
  all(table(anticlusters, schaper2019$room) == matrix(c(40, 4, 4, 40, 4, 4), ncol = 2))
)

# Test problematic cases

expect_error(
  anticlustering(
    schaper2019[, 3:6],
    K = c(2, 4),
    categories = schaper2019$room
  ), 
  pattern = "misspecified"
)

expect_error(
  anticlustering(
    schaper2019[, 3:6],
    K = c(1, 2, 3, 4, 5, 100000),
    categories = schaper2019$room
  ), 
  pattern = "misspecified"
)

# too many groups
expect_error(
  anticlustering(
    schaper2019[, 3:6],
    K = nrow(schaper2019)
  ), 
  pattern = "smaller than"
)

# preclustering not possible
expect_error(
  anticlustering(
    schaper2019[, 3:6],
    K = c(24, 24, 48),
    preclustering = TRUE
  ), 
  pattern = "preclustering"
)

