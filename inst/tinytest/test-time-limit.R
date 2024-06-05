
### Test that time limit argument works for all solvers

## Diversity

N <- 40
M <- 5
K <- 4

data <- matrix(rnorm(N*M), ncol = M)

expect_error(
  optimal_anticlustering(data, K, objective = "diversity", solver = "lpSolve", time_limit = 1),
  pattern = "time limit"
)

expect_error(
  optimal_anticlustering(data, K, objective = "diversity", solver = "glpk", time_limit = 1),
  pattern = "time limit"
)

expect_error(
  optimal_anticlustering(data, K, objective = "diversity", solver = "symphony", time_limit = 1),
  pattern = "time limit"
)

## Dispersion

N <- 400
M <- 5
K <- 8

data <- matrix(rnorm(N*M), ncol = M)

expect_error(
  optimal_anticlustering(data, K, objective = "dispersion", solver = "lpSolve", time_limit = 1),
  pattern = "time limit"
)

expect_error(
  opt <- optimal_anticlustering(data, K, objective = "dispersion", solver = "glpk", time_limit = 1),
  pattern = "time limit"
)

expect_error(
  optimal_anticlustering(data, K, objective = "dispersion", solver = "symphony", time_limit = 1),
  pattern = "time limit"
)

