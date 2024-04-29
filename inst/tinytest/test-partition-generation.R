

context("Generate all partitions")
library("anticlust")

test_that("generate_partitions function correctly removes redundant partitions", {
  for (K in 2:4) {
    for (N in 4:10) {
      if (N %% K != 0) {
        next
      }
      permutations <- generate_partitions(N, K, TRUE)
      partitions   <- generate_partitions(N, K, FALSE)
      # remove duplicates from permutations
      permutations <- lapply(permutations, order_cluster_vector)
      permutations <- permutations[!duplicated(permutations)]
      expect_equal(permutations, partitions)
    }
  }
})

test_that("analytical solution and generative functions result in same number of partitions", {
  K <- 2
  for (N in seq(4, 18, 2)) {
    partitions <- generate_partitions(N, K, FALSE)
    analytical_n <- n_partitions(N, K)
    expect_equal(analytical_n, length(partitions))
  }
  K <- 3
  for (N in seq(6, 12, 3)) {
    partitions <- generate_partitions(N, K, FALSE)
    analytical_n <- n_partitions(N, K)
    expect_equal(analytical_n, length(partitions))
  }
  K <- 4
  for (N in c(8, 12)) {
    partitions <- generate_partitions(N, K, FALSE)
    analytical_n <- n_partitions(N, K)
    expect_equal(analytical_n, length(partitions))
  }
  K <- 5
  for (N in c(5, 10)) {
    partitions <- generate_partitions(N, K, FALSE)
    analytical_n <- n_partitions(N, K)
    expect_equal(analytical_n, length(partitions))
  }
})
