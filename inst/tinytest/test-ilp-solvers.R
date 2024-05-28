
### Test solvers for maximum diversity

N <- 14
M <- 5
K <- 2

data <- matrix(rnorm(N*M), ncol = M)

lpsolve <- optimal_anticlustering(data, K, objective = "diversity", solver = "lpSolve")
val1 <- diversity_objective(data, lpsolve)

glpk <- optimal_anticlustering(data, K, objective = "diversity", solver = "glpk")
val2 <- diversity_objective(data, glpk)

symphony <- optimal_anticlustering(data, K, objective = "diversity", solver = "symphony") # symphony is slowest for max diversity!
val3 <- diversity_objective(data, symphony)

default <- optimal_anticlustering(data, K, objective = "diversity") # Test default
val4 <- diversity_objective(data, default)

val5 <- diversity_objective(data, anticlustering(data, K = K))

expect_equal(val1, val2)
expect_equal(val1, val3)
expect_equal(val1, val4)
expect_true(val1 >= val5)


### Test solvers for maximum dispersion

x <- schaper2019[, 3:6]
start <- Sys.time()
val1 <- optimal_dispersion(x, K = 3, solver = "lpSolve")$dispersion # here lpSolve is slower
Sys.time() - start

x <- schaper2019[, 3:6]
start <- Sys.time()
val2 <- optimal_dispersion(x, K = 3, solver = "glpk")$dispersion
Sys.time() - start

x <- schaper2019[, 3:6]
start <- Sys.time()
val3 <- optimal_dispersion(x, K = 3, solver = "symphony")$dispersion # here symphony is fastest
Sys.time() - start

expect_equal(val1, val2)
expect_equal(val1, val3)

