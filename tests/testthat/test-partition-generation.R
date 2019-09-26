

context("Generate all partitions")
library("anticlust")

test_that("generate_partitions function correctly removes redundant partitions", {
  for (K in 2:4) {
    for (N in 4:10) {
      if (N %% K != 0) {
        next
      }
      permutations <- generate_partitions(K, N, TRUE)
      partitions   <- generate_partitions(K, N, FALSE)
      # remove duplicates from permutations
      permutations <- lapply(permutations, order_cluster_vector)
      permutations <- permutations[!duplicated(permutations)]
      expect_equal(permutations, partitions)
    }
  }
})
