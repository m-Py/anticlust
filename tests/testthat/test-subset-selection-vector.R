


context("Initialize vector for subset selection exchange partners")
library("anticlust")

test_that("Test initialization of vector for subset selection (includes NAs)", {

  ## Test input via N/K argument
  clusters <- subset_selection(N = 12, K = 3, n = c(2, 4, 2))
  expect_equal(all(table(clusters) == c(2, 4, 2)), TRUE)
  expect_equal(length(clusters) == 12, TRUE)

  clusters <- subset_selection(N = 15, K = 3, n = c(1, 2, 3))
  expect_equal(all(table(clusters) == c(1, 2, 3)), TRUE)
  expect_equal(length(clusters) == 15, TRUE)

  clusters <- subset_selection(N = 20, K = 4, n = c(2, 3, 4, 5))
  expect_equal(all(table(clusters) == c(2, 3, 4, 5)), TRUE)
  expect_equal(length(clusters) == 20, TRUE)

  clusters <- subset_selection(N = 20, K = 4, n = c(5, 5, 5, 5))
  expect_equal(all(table(clusters) == c(5, 5, 5, 5)), TRUE)
  expect_equal(length(clusters) == 20, TRUE)

  ## Test N/K input with variyng n and no NA (different sample sizes for stimuli!)
  clusters <- subset_selection(n = c(9, 5, 3))
  expect_equal(all(table(clusters) == c(9, 5, 3)), TRUE)
  expect_equal(length(clusters) == sum(c(9, 5, 3)), TRUE)

  ## test input via `groups` argument
  ## Per iris species, select five plants. Use random order to test robustness
  groups <- sample(iris[, 5])
  clusters <- subset_selection(n = c(5), groups = groups)
  expect_equal(
    all(table(groups, clusters) == matrix(c(5, 0, 0, 0, 5, 0, 0, 0, 5), ncol = 3)),
    TRUE
  )

  ## Use vector argument for `n`
  clusters <- subset_selection(n = c(5, 10, 15), groups = groups)
  expect_equal(
    all(table(groups, clusters) == matrix(c(5, 0, 0, 0, 10, 0, 0, 0, 15), ncol = 3)),
    TRUE
  )
})
