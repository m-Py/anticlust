

context("Test computation of objective values")
library("anticlust")

test_that("objective value for variance criterion is computed correctly", {
  for (m in 1:4) {
    m_features <- m
    ## vary number of anticlusters
    for (p in 2:5) {
      p_anticlusters <- p
      n_elements <- p * 3 # n-elements must be multiplier of p
      features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
      cl <- kmeans(features, p_anticlusters)
      obj_kmeans <- cl$tot.withinss
      obj_mine  <- variance_objective(features, cl$cluster)
      expect_equal(obj_kmeans, obj_mine)
    }
  }
})

test_that("objective value for distance criterion is computed correctly", {
  conditions <- expand.grid(m = 1:4, p = 2:3)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    p_anticlusters <- conditions[k, "p"]
    n_elements <- p_anticlusters * 3 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- as.matrix(dist(features))
    ilp <- anticlustering_ilp(distances, p_anticlusters)
    solution <- solve_ilp(ilp, "min")
    anticlusters <- ilp_to_groups(solution, n_elements)
    expect_equal(solution$obj, diversity_objective_(anticlusters, features))
  }
})
