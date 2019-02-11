

context("Test computation of objective values")
library("anticlust")

test_that("output of objective value functions has correct structure", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (i in nrow(conditions)) {
    m_features <- conditions[i, "m"]
    p_anticlusters <- conditions[i, "p"]
    n_elements <- p_anticlusters * 5 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    # Precluster cases
    n_preclusters <- nrow(features) / p_anticlusters
    preclusters <- equal_sized_kmeans(features, n_preclusters)
    # Use preclustering as resticting information in anticlustering
    anticlusters <- heuristic_anticlustering(features, preclusters, nrep = 100,
                                             objective = "distance",
                                             method = "sampling")
    dist_obj <- get_objective(features, anticlusters, "distance")
    var_obj <- get_objective(features, anticlusters, "variance")
    expect_equal(mode(dist_obj), "numeric")
    expect_equal(length(dist_obj), 1)
    expect_equal(dist_obj > 0, TRUE)
    expect_equal(mode(var_obj), "numeric")
    expect_equal(length(var_obj), 1)
    expect_equal(var_obj > 0, TRUE)
  }
})


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
      obj_mine  <- obj_value_variance(features, cl$cluster)
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
    distances <- dist(features)
    solver <- "Rglpk"
    ilp <- anticlustering_ilp(distances, p_anticlusters, solver)
    solution <- solve_ilp(ilp, solver, "min")
    anticlusters <- ilp_to_groups(ilp, solution)
    expect_equal(solution$obj, get_objective(features, anticlusters, "distance"))
  }
})
