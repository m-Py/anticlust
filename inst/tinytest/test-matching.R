

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
  matches <- matching(data, p = p)
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
  n <- sample(40:100, size = 1)
  data <- matrix(rnorm(n * m), ncol = m)
  p <- sample(2:5, size = 1)
  n_groups <- sample(2:5, size = 1)
  groups <- sample(1:n_groups, size = n, replace = TRUE)
  while (any(table(groups) < p) || length(unique(groups)) != n_groups) {
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
  matches <- matching(dist(data), p = p, match_within = groups)
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

test_that("Matching behaves correctly when combining `match_within` and `match_between`", {
  # generate some random data
  p <- sample(2:5, size = 1)
  tab <- p - 1
  while (any(tab < p)) {
    m <- sample(1:4, size = 1)
    n <- sample(20:200, size = 1)
    data <- matrix(rnorm(n * m), ncol = m)
    n_groups_within <- sample(2:5, size = 1)
    groups_between <- sample(1:p, size = n, replace = TRUE)
    groups_within <- sample(1:n_groups_within, size = n, replace = TRUE)
    tab <- table(groups_between, groups_within)
  }
  
  # feature input
  # test that matches are of size p
  matches <- matching(
    data, 
    p = p, 
    match_between = groups_between,
    match_within = groups_within
  )
  expect_true(all(table(matches) == p))
  # predict number of matches
  tab <- table(groups_within, groups_between)
  # for each group_within, there are as many matches as the minimum groups_between
  n_matches <- sum(apply(tab, 1, min))
  expect_equal(n_matches, max(matches, na.rm = TRUE))
  # all not-matched elements have NA
  expect_equal(n - (n_matches * p), sum(is.na(matches)))
  # check the balancing across all grouping variables
  tab <- table(matches, groups_within, groups_between)
  expect_true(all(tab %in% c(0, 1)))
  expect_true(all(apply(tab, 3, rowSums) == 1))
  
  # feature input
  # repeat the same for distance input
  matches <- matching(
    dist(data), 
    p = p, 
    match_between = groups_between,
    match_within = groups_within
  )
  expect_true(all(table(matches) == p))
  # predict number of matches
  tab <- table(groups_within, groups_between)
  # for each group_within, there are as many matches as the minimum groups_between
  n_matches <- sum(apply(tab, 1, min))
  expect_equal(n_matches, max(matches, na.rm = TRUE))
  # all not-matched elements have NA
  expect_equal(n - (n_matches * p), sum(is.na(matches)))
  # check the balancing across all grouping variables
  tab <- table(matches, groups_within, groups_between)
  expect_true(all(tab %in% c(0, 1)))
  expect_true(all(apply(tab, 3, rowSums) == 1))
})

# this test is monster code; it tests that the first target element
# is selected as it should be and paired with its nearest neighbour,
# for different input. Fun fact: I made a lot of mistakes implementing
# this test, but the code itself was fine all the time ... 
test_that("Algorithm matches the element it should match", {
  # generate some random data
  m <- sample(1:4, size = 1)
  n <- sample(40:100, size = 1)
  data <- matrix(rnorm(n * m), ncol = m)
  p <- sample(2:5, size = 1)
  
  for (most_extreme in c(TRUE, FALSE)) {
    if (most_extreme == TRUE) {
      FUN <- which.max
    } else {
      FUN <- which.min
    }
    for (dat in c("features", "distance")) {
      if (dat == "distances") {
        df <- dist(data)
        distances <- as.matrix(data)
      } else {
        df <- data
        distances <- as.matrix(dist(data))
      }
      # vary if I use p or matches_between as input
      for (j in c("p", "groups")) {
        if (j == "p") {
          matches <- matching(df, p = p, match_extreme_first = most_extreme)
          first_target <- FUN(distances_from_centroid(data))
          target_group <- FALSE
        } else {
          groups <- to_numeric(sample(1:p, size = n, replace = TRUE))
          while (any(table(groups) < 2) || length(unique(table(groups))) == 1) {
            groups <- sample(1:p, size = n, replace = TRUE)
          }
          groups <- merge_into_one_variable(groups)
          matches <- matching(df, match_between = groups, match_extreme_first = most_extreme)
          centroid_distances <- distances_from_centroid(data)
          smallest_group <- which.min(table(groups))
          # select index of element in smallest group based on distance to centroid
          target_group <- groups == smallest_group
          centroid_distances[!target_group] <- NA
          first_target <- FUN(centroid_distances)
        }
        target_distances <- distances[, first_target]
        target_distances[target_group] <- NA # only select from «other» groups
        closest <- min(target_distances[-first_target], na.rm = TRUE)
        closest_neighbours <- which(target_distances == closest)
        # ensure that at least one of the neighbours with minimum distance
        # is in the same match
        expect_true(closest_neighbours %in% which(matches == matches[first_target]))
      }
    }
  }
})
