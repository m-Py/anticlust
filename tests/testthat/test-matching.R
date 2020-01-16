

context("Test matching function")
library("anticlust")

test_that("Matching function behaves correctly with `p` argument", {
  # generate some random data
  m <- sample(1:4, size = 1)
  n <- sample(10:100, size = 1)
  data <- matrix(rnorm(n * m), ncol = m)
  p <- sample(2:5, size = 1)
  
  # feature input
  # test that matches are of size p
  matches <- matching(data, p = p)
  expect_true(all(table(matches) == p))
  # non-fitting elements have NA
  expect_equal(sum(is.na(matches)), n %% p)
  # are all matches objective sorted by distance objective?
  objectives <- sapply(
    1:max(matches, na.rm = TRUE), 
    function(x) sum(dist(data[!is.na(matches) & matches == x, ]))
  )
  expect_true(!is.unsorted(objectives))
  
  # repeat the above for distance input
  data <- as.matrix(dist(data))
  matches <- matching(distances = data, p = p)
  expect_true(all(table(matches) == p))
  expect_equal(sum(is.na(matches)), n %% p)
  objectives <- sapply(
    1:max(matches, na.rm = TRUE), 
    function(x) sum(as.dist(data[!is.na(matches) & matches == x, !is.na(matches) & matches == x]))
  )
  expect_true(!is.unsorted(objectives))
})

test_that("Matching function behaves correctly with `match_between` argument", {
  # generate some random data
  m <- sample(1:4, size = 1)
  n <- sample(20:100, size = 1)
  data <- matrix(rnorm(n * m), ncol = m)
  p <- sample(2:5, size = 1)
  groups <- sample(1:p, size = n, replace = TRUE)
  while (any(table(groups) < 2)) {
    groups <- sample(1:p, size = n, replace = TRUE)
  }
  
  # feature input
  # test that matches are of size 
  matches <- matching(data, match_between = groups)
  expect_true(all(table(matches) == p))
  # non-fitting elements have NA
  n_matched <- (min(table(groups)) * p) # how many elements were matched
  expect_equal(sum(is.na(matches)), n - n_matched)
  # are all matches objective sorted by distance objective?
  objectives <- sapply(
    1:max(matches, na.rm = TRUE), 
    function(x) sum(dist(data[!is.na(matches) & matches == x, ]))
  )
  expect_true(!is.unsorted(objectives))
  
  # repeat the above for distance input
  data <- as.matrix(dist(data))
  matches <- matching(data, match_between = groups)
  expect_true(all(table(matches) == p))
  # non-fitting elements have NA
  n_matched <- (min(table(groups)) * p) # how many elements were matched
  expect_equal(sum(is.na(matches)), n - n_matched)
  objectives <- sapply(
    1:max(matches, na.rm = TRUE), 
    function(x) sum(as.dist(data[!is.na(matches) & matches == x, !is.na(matches) & matches == x]))
  )
  expect_true(!is.unsorted(objectives))
})


test_that("Matching function behaves correctly with `match_within` argument", {
  # generate some random data
  m <- sample(1:4, size = 1)
  n <- sample(20:100, size = 1)
  data <- matrix(rnorm(n * m), ncol = m)
  p <- sample(2:5, size = 1)
  n_groups <- sample(2:5, size = 1)
  groups <- sample(1:n_groups, size = n, replace = TRUE)
  while (any(table(groups) < p)) {
    groups <- sample(1:p, size = n, replace = TRUE)
  }
  
  # feature input
  # test that matches are of size p
  matches <- matching(data, p = p, match_within = groups)
  expect_true(all(table(matches) == p))
  # test that all matches are within a category
  tab <- table(matches, groups)
  expect_true(all(apply(tab, 1, function(x) sum(x == p)) == 1))
  expect_true(all(apply(tab, 1, function(x) sum(x == 0)) == (n_groups - 1)))
  # non-fitting elements have NA
  not_matched <- sum(table(groups) %% p) # how many elements were matched
  expect_equal(sum(is.na(matches)), not_matched)
  # are all matches objective sorted by distance objective?
  objectives <- sapply(
    1:max(matches, na.rm = TRUE), 
    function(x) sum(dist(data[!is.na(matches) & matches == x, ]))
  )
  expect_true(!is.unsorted(objectives))
  
  # repeat the above for distance input
  # test that matches are of size p
  matches <- matching(distances = dist(data), p = p, match_within = groups)
  expect_true(all(table(matches) == p))
  tab <- table(matches, groups)
  expect_true(all(apply(tab, 1, function(x) sum(x == p)) == 1))
  expect_true(all(apply(tab, 1, function(x) sum(x == 0)) == (n_groups - 1)))
  not_matched <- sum(table(groups) %% p) 
  expect_equal(sum(is.na(matches)), not_matched)
  objectives <- sapply(
    1:max(matches, na.rm = TRUE), 
    function(x) sum(dist(data[!is.na(matches) & matches == x, ]))
  )
  expect_true(!is.unsorted(objectives))
})
