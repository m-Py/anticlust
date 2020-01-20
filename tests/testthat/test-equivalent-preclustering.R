
context("Differenct preclustering interfaces should produce the same outputs")
library("anticlust")

test_that("different options for preclustering have the same result - variance objective", {
  for (M in 1:4) {
    for (K in 2:5) {
      N <- K * 3
      features <- matrix(rnorm(N * M), ncol = M)

      # determine a random seed to make the exchange method reproducible
      seed <- sample(1000, size = 1)
      set.seed(seed)

      ## First option
      # Set `preclustering = TRUE`
      ac1 <- anticlustering(
        features,
        K = K,
        objective = "variance",
        preclustering = TRUE
      )

      ## Second option
      # Call `balanced_clustering` and use output as `categories` argument
      preclusters <- balanced_clustering(
        features,
        K = N / K
      )
      set.seed(seed)
      ac2 <- anticlustering(
        features,
        K = K,
        objective = "variance",
        categories = preclusters
      )

      ## Third option
      # Use `fast_anticlustering` function
      set.seed(seed)
      ac3 <- fast_anticlustering(
        features,
        K = K,
        categories = preclusters
      )

      expect_equal(all(ac1 == ac2), TRUE)
      expect_equal(all(ac2 == ac3), TRUE)
    }
  }
})


test_that("different options for preclustering have the same result - distance objective", {
  for (M in 1:4) {
    for (K in 2:5) {
      N <- K * 3
      features <- matrix(rnorm(N * M), ncol = M)

      # determine a random seed to make the exchange method reproducible
      seed <- sample(1000, size = 1)
      set.seed(seed)

      ## First option
      # Set `preclustering = TRUE`
      ac1 <- anticlustering(
        features,
        K = K,
        objective = "distance",
        preclustering = TRUE
      )

      ## Second option
      # Call `balanced_clustering` and use output as `categories` argument
      preclusters <- balanced_clustering(
        features,
        K = N / K
      )
      set.seed(seed)
      ac2 <- anticlustering(
        features,
        K = K,
        objective = "distance",
        categories = preclusters
      )

      ## Third option
      # Use distance input (and categories for preclusters because the
      # preclustering algorithm is slightly different for distance input)
      set.seed(seed)
      ac3 <- anticlustering(
        dist(features),
        K = K,
        objective = "distance",
        categories = preclusters
      )

      expect_equal(all(ac1 == ac2), TRUE)
      expect_equal(all(ac2 == ac3), TRUE)
    }
  }
})
