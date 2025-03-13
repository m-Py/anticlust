
### Test solvers for maximum diversity
library(anticlust)
library(tinytest)
N <- 12
M <- 5
K <- 2

data <- matrix(rnorm(N*M), ncol = M)
start <- Sys.time()
lpsolve <- optimal_anticlustering(data, K, objective = "diversity", solver = "lpSolve")
Sys.time() - start
val1 <- diversity_objective(data, lpsolve)

start <- Sys.time()
glpk <- optimal_anticlustering(data, K, objective = "diversity", solver = "glpk")
Sys.time() - start
val2 <- diversity_objective(data, glpk)

start <- Sys.time()
symphony <- optimal_anticlustering(data, K, objective = "diversity", solver = "symphony") # symphony is slowest for max diversity!
Sys.time() - start
val3 <- diversity_objective(data, symphony)

expect_equal(val1, val2)
expect_equal(val1, val3)

if (requireNamespace("gurobi", quietly = TRUE)) {
  gurobi <- optimal_anticlustering(data, K, objective = "diversity", solver = "gurobi") 
  val4 <- diversity_objective(data, gurobi)
  expect_equal(val1, val4)
}

### Test solvers for maximum dispersion

x <- schaper2019[, 3:6]
start <- Sys.time()
val1 <- optimal_dispersion(x, K = 3, solver = "lpSolve")$dispersion # here lpSolve is slower
Sys.time() - start
start <- Sys.time()
val2 <- optimal_dispersion(x, K = 3, solver = "glpk")$dispersion
Sys.time() - start
start <- Sys.time()
val3 <- optimal_dispersion(x, K = 3, solver = "symphony")$dispersion
Sys.time() - start
expect_equal(val1, val2)
expect_equal(val1, val3)

if (requireNamespace("gurobi", quietly = TRUE)) {
  val4 <- optimal_dispersion(x, K = 3, solver = "gurobi")$dispersion
  expect_equal(val1, val4)
}



### Test solvers for balanced clustering (i.e., reversed maximum - minimum - diversity)

N <- 20
M <- 5
K <- 4

data <- matrix(rnorm(N*M), ncol = M)

lpsolve <- balanced_clustering(data, K, method = "ilp", solver = "lpSolve")
val1 <- diversity_objective(data, lpsolve)

glpk <- balanced_clustering(data, K, method = "ilp", solver = "glpk")
val2 <- diversity_objective(data, glpk)

symphony <- balanced_clustering(data, K, method = "ilp", solver = "symphony")
val3 <- diversity_objective(data, symphony)

expect_equal(val1, val2)
expect_equal(val1, val3)

if (requireNamespace("gurobi", quietly = TRUE)) {
  gurobi <- balanced_clustering(data, K, method = "ilp", solver = "gurobi")
  val4 <- diversity_objective(data, gurobi)
  expect_equal(val1, val4)
}

### Test solvers for weighted cluster editing (i.e., reversed max diversity without group restrictions)

features <- swiss[sample(nrow(swiss), size = 20), ]
features <- swiss
distances <- dist(scale(features))
agreements <- ifelse(as.matrix(distances) < median(distances), 1, -1)

start <- Sys.time()
lpsolve <- wce(agreements, solver = "lpSolve")
Sys.time() - start
val1 <- diversity_objective(agreements, lpsolve)

start <- Sys.time()
glpk <- wce(agreements, solver = "glpk")
Sys.time() - start
val2 <- diversity_objective(agreements, glpk)

start <- Sys.time()
symphony <- wce(agreements, solver = "symphony")
Sys.time() - start
val3 <- diversity_objective(agreements, symphony)

expect_equal(val1, val2)
expect_equal(val1, val3)

if (requireNamespace("gurobi", quietly = TRUE)) {
  gurobi <- wce(agreements, solver = "gurobi")
  val4 <- diversity_objective(agreements, gurobi)
  expect_equal(val1, val4)
}

