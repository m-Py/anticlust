library(anticlust)
library(tinytest)

set.seed(123)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)

result_cluster1 <- three_phase_search_anticlustering(distances, K, N)
result_cluster2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
result_cluster3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve")

diversity1 <- diversity_objective(distances, result_cluster1)
diversity2 <- diversity_objective(distances, result_cluster2)
diversity3 <- diversity_objective(distances, result_cluster3)

expect_true(diversity3 >= diversity1)
expect_true(diversity3 >= diversity2)

### Test more clusters ###

N <- 12
M <- 2
K <- 3

dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)

result_cluster1 <- anticlust:::three_phase_search_anticlustering(distances, K, N)
result_cluster2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
result_cluster3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve", time_limit = 20)

diversity1 <- diversity_objective(distances, result_cluster1)
diversity2 <- diversity_objective(distances, result_cluster2)
diversity3 <- diversity_objective(distances, result_cluster3)

expect_true(diversity3 >= diversity1)
expect_true(diversity3 >= diversity2)


## test cluster vector with different sizes

N2 <- 11
M2 <- 2
K2 <- 3
clusters <- c(3,3,5)
dat2 <- matrix(rnorm(N2 * M2), ncol = M2)
distances2 <- dist(dat2)

result_cluster <- anticlust:::three_phase_search_anticlustering(distances2, K2, N2, clusters=clusters)
table_cluseters <- table(result_cluster)
table_cluseters
expect_true(all(table_cluseters == clusters))

## test cluster vector swith higher cluster size

N2 <- 140
M2 <- 2
K2 <- 4
clusters <- c(20,40,30,50)
dat2 <- matrix(rnorm(N2 * M2), ncol = M2)
distances2 <- dist(dat2)

result_cluster <- anticlust:::three_phase_search_anticlustering(distances2, K2, N2, clusters=clusters)
table_cluseters <- table(result_cluster)
table_cluseters
expect_true(all(table_cluseters == clusters))


## test cluster vector throughs error

N2 <- 11
M2 <- 2
K2 <-2
clusters <- c(3,3,4)
dat2 <- matrix(rnorm(N2 * M2), ncol = M2)
distances2 <- dist(dat2)

expect_error(anticlust:::three_phase_search_anticlustering(distances2, K2, N2, clusters=clusters))


# Test problematic cases

expect_error(
  anticlust:::three_phase_search_anticlustering(distances2, K2, N2, clusters=c(1)),
  pattern = "len"
)

expect_error(
  anticlust:::three_phase_search_anticlustering(distances2, K2, N2, alpha=-1),
  pattern = "alpha must be greater than 0"
)

expect_error(
  anticlust:::three_phase_search_anticlustering(distances2, K2, N2, eta_max=1.2),
  pattern = "must be integer"
)

expect_error(
  anticlust:::three_phase_search_anticlustering(distances2, K2, N2, objective = "something"),
  pattern = "Argument objective can either be set to 'diversity' or 'dispersion"
)


### Test dispersion ###


set.seed(123)

N <- 12
M <- 2
K <- 3

dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)

result_cluster1 <- anticlust:::three_phase_search_anticlustering(dat, K, N, objective = "dispersion", number_iterations = 1, beta_max = 5)

result_cluster2 <- optimal_anticlustering(distances, objective = "dispersion", K=K, solver = "lpSolve")
results <- optimal_dispersion(dat, K=K, solver = "symphony")$dispersion

dispersion1 <- dispersion_objective(distances, result_cluster1)
dispersion2 <- dispersion_objective(distances, result_cluster2)

