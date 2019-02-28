library("anticlust")

context("Simulated annealing")

test_that("Candidate generation for simulated annealing works correctly", {
  for (i in 2:4) {
    original <- rep(1:i, 10)
    for (j in c(TRUE, FALSE)) {
      next_ <- next_candidate(original, j)
      ## Exactly two elements have to differ
      expect_equal(sum(original != next_), 2)
      ## Still legal number of elements per cluster
      expect_equal(all(table(next_) == table(next_)[1]), TRUE)
    }
  }
})
