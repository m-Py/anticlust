
context("C implementation")
library("anticlust")

test_that("C implemenation has same output as R implementation", {
  set.seed(123)
  N <- sample(40:80, size = 1)
  M <- sample(1:5, size = 1)
  K <- sample(2:5, size = 1)
  features <- matrix(rnorm(N * M), ncol = M)
  clusters <- initialize_clusters(N, K, NULL)
  cl1 <- anticlustering(features, clusters, objective = "variance")
  cl2 <- anticlustering(features, clusters, variance_objective_)
  expect_true(all(cl1 == cl2))

  # now also use a categorical constraint
  set.seed(71207)
  C <- sample(2:5, size = 1)
  categories <- sample(C, size = N, replace = TRUE)
  clusters <- initialize_clusters(N, K, categories)
  cl1 <- anticlustering(
    features, 
    clusters, 
    objective = "variance", 
    categories = categories
  )
  cl2 <- anticlustering(
    features, 
    clusters, 
    objective = variance_objective_, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))

  # TODO: include a test for the edge cases (i.e., categories with 1 or N members)
  
  # Use diversity objective
  set.seed(123)
  clusters <- initialize_clusters(N, K, NULL)
  cl1 <- anticlustering(features, clusters, objective = "distance")
  cl2 <- anticlustering(features, clusters, distance_objective_)
  expect_true(all(cl1 == cl2))
  
  # now also use a categorical constraint
  set.seed(71207)
  C <- sample(2:5, size = 1)
  categories <- sample(C, size = N, replace = TRUE)
  clusters <- initialize_clusters(N, K, categories)
  cl1 <- anticlustering(
    features, 
    clusters, 
    objective = "distance", 
    categories = categories
  )
  cl2 <- anticlustering(
    features, 
    clusters, 
    objective = distance_objective_, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))

})
