
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

