

library("anticlust")

# k-means and diversity are equivalent for equal sized groups
N <- 50
M <- 5
K <- 5
set.seed(123)
data <- matrix(rnorm(N * M), ncol = M)
init <- anticlust:::initialize_clusters(N, K, NULL)

# k-means
g1 <- anticlustering(data, K = init, objective = "variance")
g2 <- anticlustering(dist(data)^2, K = init)
# this can go wrong due to floating point errors, therefore I use a seed
expect_true(all(g1 == g2))

# k-plus 
g1 <- anticlustering(data, K = init, objective = "kplus")
g2 <- anticlustering(dist(kplus_moment_variables(data, 2, FALSE))^2, K = init)
expect_true(all(g1 == g2))

# unequal-sized groups
# k-means
K <- c(5, 10, 15, 20)
init <- anticlust:::initialize_clusters(N, K, NULL)
g1 <- anticlustering(data, K = init, objective = "variance")
g2 <- anticlustering(dist(data)^2, K = init)
expect_false(all(g1 == g2))

# needs "average" diversity
g3 <- anticlustering(dist(data)^2, K = init, objective = "average-diversity")
expect_true(all(g1 == g3))

# k-plus 
g1 <- anticlustering(data, K = init, objective = "kplus")
g2 <- anticlustering(dist(kplus_moment_variables(data, 2, FALSE))^2, K = init)
expect_false(all(g1 == g2))

# needs "average" diversity
g3 <- anticlustering(dist(kplus_moment_variables(data, 2, FALSE))^2, K = init, objective = "average-diversity")
expect_true(all(g1 == g3))
