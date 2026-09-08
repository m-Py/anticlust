library(anticlust)
library(tinytest)

set.seed(123)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)

result_cluster1 <- feasible_and_infeasible_region_search_anticlustering(distances, K, N)
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

result_cluster1 <- feasible_and_infeasible_region_search_anticlustering(distances, K, N)
result_cluster2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
result_cluster3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve", time_limit = 20)

diversity1 <- diversity_objective(distances, result_cluster1)
diversity2 <- diversity_objective(distances, result_cluster2)
diversity3 <- diversity_objective(distances, result_cluster3)

expect_true(diversity3 >= diversity1)
expect_true(diversity3 >= diversity2)


## test cluster vector swith higher cluster size

N2 <- 140
M2 <- 2
K2 <- 4
clusters <- c(20,40,30,50)
dat2 <- matrix(rnorm(N2 * M2), ncol = M2)
distances2 <- dist(dat2)

result_cluster <- feasible_and_infeasible_region_search_anticlustering(distances2, K2, N2, clusters=clusters)
table_clusters <- table(result_cluster)
table_clusters
expect_true(all(table_clusters == clusters))



# Test problematic cases

expect_error(
  feasible_and_infeasible_region_search_anticlustering(distances2, K2, N2, clusters=c(1)),
  pattern = "len"
)

## Larger data set

N <- 200
M <- 5
K <- 10

data <- matrix(rnorm(N*M), ncol = M)

g1 <- anticlustering(data, K = K)
g2 <- anticlustering(data, K = K, method = "local-maximum")
g3 <- feasible_and_infeasible_region_search_anticlustering(data, K, N)
g4 <- three_phase_search_anticlustering(data, K, N)

diversity_objective(data, g1)
diversity_objective(data, g2)
diversity_objective(data, g3)
diversity_objective(data, g4)


## Also test this algorithm via interface in anticlustering()

data <- brunel2025[, -(1:2)]

a1 <- anticlustering(data, K = 13, standardize = TRUE, objective = "kplus")
a2 <- anticlustering(data, K = 13, method = "fifr", standardize = TRUE, objective = "kplus")
a3 <- anticlustering(data, K = 13, method = "3phase", standardize = TRUE, objective = "kplus")
mean_sd_tab(brunel2025[, -(1:4)], a1, return_diff = TRUE)
mean_sd_tab(brunel2025[, -(1:4)], a2, return_diff = TRUE)
mean_sd_tab(brunel2025[, -(1:4)], a3, return_diff = TRUE)

table(data$target_word_emotionality, a1)
table(data$target_word_emotionality, a2)
table(data$target_word_emotionality, a3)

# compare objectives, need actual features that were used in anticlustering:
ff <- anticlust:::get_anticlustering_features(data, "kplus", TRUE)
variance_objective(ff, a1)
variance_objective(ff, a2)
variance_objective(ff, a3)

