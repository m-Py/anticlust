library(tinytest)
library(anticlust)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)

# just make sure that all runs okay
optimal_anticlustering(dat, K, "diversity", solver = "glpk")
optimal_anticlustering(dat, K, "kplus", solver = "glpk")
optimal_anticlustering(dat, K, "variance", solver = "glpk")
optimal_anticlustering(dat, K, "dispersion", solver = "glpk") 

expect_error(
  optimal_anticlustering(dat, K = 5, "diversity", solver = "glpk"),
  pattern = "equal-sized groups"
)

expect_error(
  optimal_anticlustering(dat, K = c(2, 2, 2, 6), "diversity", solver = "glpk"),
  pattern = "equal-sized groups"
)

# it should work for dispersion however: 
optimal_anticlustering(dat, K = 5, "dispersion", solver = "glpk")
optimal_anticlustering(dat, K = c(2, 2, 2, 6), "dispersion", solver = "glpk")

expect_error(
  optimal_anticlustering(dat, K = K, "diversity", solver = "FOO"),
  pattern = "Argument solver"
)

expect_error(
  optimal_anticlustering(dist(dat), K = K, "variance", solver = "glpk"),
  pattern = "objective 'variance'"
)

expect_error(
  optimal_anticlustering(dist(dat), K = K, "kplus", solver = "glpk"),
  pattern = "objective 'variance'"
)


