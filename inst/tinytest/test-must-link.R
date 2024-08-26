

library(anticlust)
library(tinytest)
#set.seed(1234)
# Generate some example data, N = 128, M = 5
K <- 20
group_sizes <- 12
N <- group_sizes * K # here 128
M <- 5
data <- matrix(rnorm(N * M), ncol = M)
distances <- as.matrix(dist(data))

# Generate random patient IDs
must_link <- sample(N, replace = TRUE)

tt <- anticlustering(distances, K, must_link = must_link, method = "exchange")

# validate that must-link constraints were preserved
same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
for (i in same) {
  expect_true(all(tt[must_link == i] == tt[must_link == i][1]))
}

tt <- anticlustering(distances, K, must_link = must_link, method = "exchange", repetitions = 2)
# validate that must-link constraints were preserved
same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
for (i in same) {
  expect_true(all(tt[must_link == i] == tt[must_link == i][1]))
}

tt <- anticlustering(distances, K, must_link = must_link, method = "local-maximum")
# validate that must-link constraints were preserved
same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
for (i in same) {
  expect_true(all(tt[must_link == i] == tt[must_link == i][1]))
}

tt <- anticlustering(distances, K, must_link = must_link, method = "local-maximum", repetitions = 2)
# validate that must-link constraints were preserved
same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
for (i in same) {
  expect_true(all(tt[must_link == i] == tt[must_link == i][1]))
}

## Expect errors

expect_error(
  anticlustering(distances, K = K, must_link = must_link, categories = sample(LETTERS, size = N, replace = TRUE)),
  pattern = "categories"
)
expect_error(
  anticlustering(distances, K = K, must_link = must_link, preclustering = TRUE),
  pattern = "preclustering"
)

expect_error(
  anticlustering(data, K = K, must_link = "A"),
  pattern = "numeric"
)

expect_error(
  anticlustering(data, K = K, must_link = must_link, objective = "kplus"),
  pattern = "diversity"
)
expect_error(
  anticlustering(data, K = K, must_link = must_link, objective = "variance"),
  pattern = "diversity"
)
expect_error(
  anticlustering(data, K = K, must_link = must_link, objective = "dispersion"),
  pattern = "diversity"
)


expect_error(
  anticlustering(data, K = K, must_link = matrix(NA))
)

