

context("Exchange method updating objective")
library("anticlust")

test_that("computing distance objectives in specialized exchange method is equal to generic method", {
  for (M in 1:4) {
    for (K in 2:5) {
      for (i in 3:8) {
        N <- K * i
        features <- matrix(rnorm(N * M), ncol = M)
        distances <- as.matrix(dist(features))
        clusters <- sample(rep(1:K, N/K))

        objective <- diversity_objective_(clusters, distances)
        objective2 <- diversity_objective_(clusters, features)
        expect_equal(objective, objective2)

        # Swap two items and check objective
        to_swap <- sample(1:K, size = 2)
        swap1 <- which(clusters == to_swap[1])[1]
        swap2 <- which(clusters == to_swap[2])[1]

        # Method 1: As in specialized exchange method for anticluster editing
        selected <- selection_matrix_from_clusters(clusters)
        obj1 <- update_objective_distance(distances, selected, swap1, swap2, objective)

        # Method 2: As in generic exchange method
        obj2 <- update_objective_generic(distances, clusters, swap1, swap2, diversity_objective)
        expect_equal(obj1, obj2)
        # note, this DOES NOT always work:
        # expect_true(obj1 == obj2) # floating points =/
      }
    }
  }
})
