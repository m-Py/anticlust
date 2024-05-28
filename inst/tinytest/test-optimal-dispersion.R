
library("anticlust")

x <- schaper2019[, 3:6]
val1 <- optimal_dispersion(x, K = 3)$dispersion
val2 <- dispersion_objective(x, anticlustering(x, K = 3, method = "brusco", objective = "dispersion"))
expect_true(val1 >= val2) 


# Output format is as expected for optimal_dispersion()
N <- 30
M <- 5
K <- 3
data <- matrix(rnorm(N*M), ncol = M)
distances <- dist(data)

opt <- optimal_dispersion(distances, K = K)

groups_heuristic <- anticlustering(
  distances, 
  K = K,
  method = "local-maximum", 
  objective = "dispersion", 
  repetitions = 10
)
expect_true(dispersion_objective(distances, opt$groups) >= dispersion_objective(distances, groups_heuristic))

expect_true(all(sort(table(opt$groups)) == sort(table(anticlust:::initialize_clusters(N, K, NULL)))))


# test for unequal group sizes (groups differ by 1 at most)
N <- 28
M <- 5
K <- 3
data <- matrix(rnorm(N*M), ncol = M)
distances <- dist(data)

opt <- optimal_dispersion(distances, K = K)

groups_heuristic <- anticlustering(
  distances, 
  K = K,
  method = "local-maximum", 
  objective = "dispersion", 
  repetitions = 10
)
expect_true(dispersion_objective(distances, opt$groups) >= dispersion_objective(distances, groups_heuristic))
expect_true(all(sort(table(opt$groups)) == sort(table(anticlust:::initialize_clusters(N, K, NULL)))))


# now "truly" unequal group sizes 
K <- c(10, 12, 6)
opt <- optimal_dispersion(distances, K = K)

groups_heuristic <- anticlustering(
  distances, 
  K = K,
  method = "local-maximum", 
  objective = "dispersion", 
  repetitions = 10
)

expect_true(dispersion_objective(distances, opt$groups) >= dispersion_objective(distances, groups_heuristic))
expect_true(all(sort(table(opt$groups)) == sort(K)))

# another grouping, more extreme case
K <- c(2, 2, 24)
opt <- optimal_dispersion(distances, K = K)

groups_heuristic <- anticlustering(
  distances, 
  K = K,
  method = "local-maximum", 
  objective = "dispersion", 
  repetitions = 10
)

expect_true(dispersion_objective(distances, opt$groups) >= dispersion_objective(distances, groups_heuristic))
expect_true(all(sort(table(opt$groups)) == sort(K)))

