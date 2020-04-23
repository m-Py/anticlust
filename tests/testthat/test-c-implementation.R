
context("C implementation")
library("anticlust")

test_that("C implemenation has same output as R implementation", {
  N <- sample(40:80, size = 1)
  M <- sample(1:5, size = 1)
  K <- sample(2:5, size = 1)
  features <- matrix(rnorm(N * M), ncol = M)
  clusters <- initialize_clusters(N, K, NULL)
  cl1 <- fanticlust(features, clusters)
  cl2 <- anticlustering(features, clusters, objective = "variance")
  expect_true(all(cl1 == cl2))
  
  # now also use a categorical constraint
  C <- sample(2:5, size = 1)
  categories <- sample(C, size = N, replace = TRUE)
  clusters <- initialize_clusters(N, K, categories)
  cl1 <- fanticlust(features, clusters, categories)
  cl2 <- anticlustering(
    features, 
    clusters, 
    objective = "variance", 
    categories = categories
  )
  expect_true(all(cl1 == cl2))
  
  # TODO: include a test for the edge cases (i.e., categories with 1 or N members)
  
})
