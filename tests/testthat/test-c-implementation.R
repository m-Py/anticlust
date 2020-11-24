
context("C implementation")
library("anticlust")

test_that("C implemenation of Anticlustering has same output as R implementation", {

  set.seed(123)
  N <- sample(40:80, size = 1)
  M <- sample(1:5, size = 1)
  K <- sample(2:5, size = 1)
  features <- matrix(rnorm(N * M), ncol = M)
  clusters <- initialize_clusters(N, K, NULL)
  cl1 <- anticlustering(features, clusters, objective = "variance")
  cl2 <- anticlustering(features, clusters, variance_objective)
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
    objective = variance_objective, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))

  # Include a test for the edge cases (i.e., categories with 1 or N members)
  categories[1] <- C + 1 
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
    objective = variance_objective, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))

  # Include a test for the edge cases (i.e., categories with 1 or N members)
  categories <- rep(1, N)
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
    objective = variance_objective, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))

  # Use diversity objective
  set.seed(123)
  clusters <- initialize_clusters(N, K, NULL)
  cl1 <- anticlustering(features, clusters, objective = "diversity")
  cl2 <- anticlustering(features, clusters, diversity_objective)
  expect_true(all(cl1 == cl2))
  
  # now also use a categorical constraint
  set.seed(71207)
  C <- sample(2:5, size = 1)
  categories <- sample(C, size = N, replace = TRUE)
  clusters <- initialize_clusters(N, K, categories)
  cl1 <- anticlustering(
    features, 
    clusters, 
    objective = "diversity", 
    categories = categories
  )
  cl2 <- anticlustering(
    features, 
    clusters, 
    objective = diversity_objective, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))
  
  # now also use preclustering
  set.seed(71207)
  C <- sample(2:5, size = 1)
  categories <- sample(C, size = N, replace = TRUE)
  clusters <- initialize_clusters(N, K, categories)
  cl1 <- anticlustering(
    features, 
    clusters, 
    objective = "diversity", 
    categories = categories
  )
  cl2 <- anticlustering(
    features, 
    clusters, 
    objective = diversity_objective, 
    categories = categories
  )
  expect_true(all(cl1 == cl2))
  
  # test dispersion objective
  set.seed(123)
  N <- 50 # number of elements
  M <- 2  # number of variables per element
  K <- 2  # number of clusters
  random_data <- matrix(rnorm(N * M), ncol = M)
  random_clusters <- sample(rep_len(1:K, N))

  # Maximize the dispersion 
  optimized_clusters <- anticlustering(
    random_data,
    K = random_clusters, 
    objective = dispersion_objective
  )
  dispersion_objective(random_data, optimized_clusters)

  optimized_clusters2 <- anticlustering(
    random_data,
    K = random_clusters, 
    objective = "dispersion"
  )
  expect_true(all(optimized_clusters == optimized_clusters2))

})
