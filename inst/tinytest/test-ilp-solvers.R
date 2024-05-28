
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
val1 <- optimal_dispersion(x, K = 3, solver = "lpSolve")$dispersion # here lpSolve is slower
val2 <- optimal_dispersion(x, K = 3, solver = "glpk")$dispersion
val3 <- optimal_dispersion(x, K = 3, solver = "symphony")$dispersion

expect_equal(val1, val2)
expect_equal(val1, val3)

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

default <- balanced_clustering(data, K, method = "ilp") # Test default
val4 <- diversity_objective(data, default)

val5 <- diversity_objective(data, balanced_clustering(data, K = K))

expect_equal(val1, val2)
expect_equal(val1, val3)
expect_equal(val1, val4)
expect_true(val1 <= val5)


### Test solvers for weighted cluster editing (i.e., reversed max diversity without group restrictions)

features <- swiss[sample(nrow(swiss), size = 20), ]

distances <- dist(scale(features))
agreements <- ifelse(as.matrix(distances) < median(distances), 1, -1)

lpsolve <- wce(agreements, solver = "lpSolve")
val1 <- diversity_objective(agreements, lpsolve)

glpk <- wce(agreements, solver = "glpk")
val2 <- diversity_objective(agreements, glpk)

symphony <- wce(agreements, solver = "symphony") # again, Symphony seems to be fastest
val3 <- diversity_objective(agreements, symphony)

expect_equal(val1, val2)
expect_equal(val1, val3)

