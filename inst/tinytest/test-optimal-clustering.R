library(tinytest)
library(anticlust)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)

# just make sure that all runs okay
balanced_clustering(dat, K, method = "ilp", solver = "glpk")

expect_error(
  balanced_clustering(dat, K = 5, method = "ilp", solver = "glpk")
)

expect_error(
  balanced_clustering(dat, K = c(2, 2, 2, 6), method = "ilp", solver = "glpk")
)

expect_error(
  balanced_clustering(dat, K = 5, method = "ilp", solver = "FOO")
)
