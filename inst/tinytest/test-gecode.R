
# test gecode solver

library(anticlust)
library(tinytest)

if (requireNamespace("gkc.gecode", quietly = TRUE)) {
  N <- 30
  M <- 6
  K <- 3

  data <- matrix(rnorm(N*M), ncol = M)
  
  start <- Sys.time()
  val1 <- optimal_dispersion(data, K = K, solver = "gecode")$dispersion
  Sys.time() - start
  
  start <- Sys.time()
  val2 <- optimal_dispersion(data, K = K, solver = "symphony")$dispersion
  Sys.time() - start
  expect_equal(val1, val2)
  
  start <- Sys.time()
  groups <- optimal_anticlustering(data, K = K, solver = "gecode", objective = "dispersion")
  Sys.time() - start
  val3 <- dispersion_objective(data, groups)
  expect_equal(val1, val3)
  
  # Use cannot_link constraints: Element 1 must not be linked with elements 2 to 10:
  cl_matrix <- matrix(c(rep(1, 9), 2:10), ncol = 2)
  cl <- anticlustering(
    schaper2019[, 3:6],
    K = 3,
    cannot_link = cl_matrix
  )
  all(cl[1] != cl[2:10])

}
