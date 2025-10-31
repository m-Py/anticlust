
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

# Function that test if constraints are valid
must_link_constraints_valid <- function(cl, must_link) {
  same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
  all_good <- c()
  for (i in same) {
    all_good <- c(all_good, all(cl[must_link == i] == cl[must_link == i][1]))
  }
  all(all_good)
}

# Generate random patient IDs
must_link <- sample(N, replace = TRUE)

# Fulfilling the must-link constraints may fail, even though the data generating 
# process here makes it very unlikely. Still we should wrap the code in a tryCatch
# block to make sure. 

tt <- tryCatch(
  anticlustering(distances, K, must_link = must_link, method = "exchange"),
  error = function(e) e
)

if (!"simpleError" %in% class(tt)) {
  expect_true(must_link_constraints_valid(tt, must_link))
}

tt <- tryCatch(
  anticlustering(distances, K, must_link = must_link, method = "exchange", repetitions = 2),
  error = function(e) e
)

if (!"simpleError" %in% class(tt)) {
  expect_true(must_link_constraints_valid(tt, must_link))
}

tt <- tryCatch(
  anticlustering(distances, K, must_link = must_link, method = "local-maximum"),
  error = function(e) e
)

if (!"simpleError" %in% class(tt)) {
  expect_true(must_link_constraints_valid(tt, must_link))
}


tt <- tryCatch(
  anticlustering(distances, K, must_link = must_link, method = "local-maximum", repetitions = 2),
  error = function(e) e
)

if (!"simpleError" %in% class(tt)) {
  expect_true(must_link_constraints_valid(tt, must_link))
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
  pattern = "length"
)
expect_error(
  anticlustering(data, K = K, must_link = must_link, objective = "dispersion"),
  pattern = "diversity"
)

expect_error(
  anticlustering(data, K = K, must_link = must_link, cannot_link = rbind(1:2, 2:3)),
  pattern = "cannot-link"
)
expect_error(
  anticlustering(data, K = K, must_link = matrix(NA))
)

# no error for unequal-sized groups with diversity
tryCatch(anticlustering(1:100, K = c(10, 10, 80), must_link = sample(100, replace = TRUE)), error = function(e) e)
# no error for unequal-sized groups with k-means/k-plus
expect_error(
  anticlustering(1:100, K = c(10, 10, 80), must_link = sample(100, replace = TRUE), objective = "variance"),
  pattern = "equal-sized"
)
expect_error(
  anticlustering(1:100, K = c(10, 10, 80), must_link = sample(100, replace = TRUE), objective = "kplus"),
  pattern = "equal-sized"
)


#try k-plus objective with schaper data set

data <- schaper2019[, 3:6]
must_link <- rep(NA, nrow(data))
must_link[1:5] <- 1
gr <- anticlustering(data, K = 4, objective = "kplus", repetitions = 100, must_link = must_link)
expect_true(all(gr[1:5] == gr[1]))
mean_sd_tab(data, gr) #nice!
