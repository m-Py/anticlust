
context("Exchange partner generation - exported function")
library("anticlust")

test_that("function returns correct results for input argument N", {
  N <- 12
  p <- 3
  partners <- generate_exchange_partners(N, p)
  expect_true(all(table(partners) == p + 1))

  # try imbalanced data
  N <- 11
  p <- 3
  partners <- generate_exchange_partners(N, p)
  expect_true(all(table(partners) >= p))
})

test_that("function returns correct results for input argument N", {
  categories <- schaper2019$room
  p <- 3
  partners <- generate_exchange_partners(categories = categories, p = p)
  tab <- table(partners)
  expect_true(all(tab == p + 1))

  # test that all is balanced across categories
  tab <- table(partners, categories)
  not_zero <- tab[tab != 0][1]
  expect_true(all(tab %in% c(0, not_zero)))
  expect_true(diff(colSums(tab)) == 0)

  # try imbalanced data
  p <- 8
  partners <- generate_exchange_partners(categories = categories, p = p)
  tab <- table(partners)
  expect_true(all(tab %in% c(p + 2, p + 1)))
})

test_that("function returns correct results for input argument similar = TRUE", {
  p <- 3 # corresponds to clusters of size 4!
  partners <- generate_exchange_partners(
    features = schaper2019[, 3:6],
    p = p,
    similar = TRUE
  )
  expect_true(all(table(partners) == p + 1))
})

test_that("function returns correct results for input argument similar = TRUE and categories", {
  p <- 3 # corresponds to clusters of size 4!
  categories <- schaper2019$room
  partners <- generate_exchange_partners(
    features = schaper2019[, 3:4],
    p = p,
    similar = TRUE,
    categories = categories
  )
  expect_true(all(table(partners) == p + 1))
  tab <- table(partners, categories)
  not_zero <- tab[tab != 0][1]
  expect_true(all(tab %in% c(0, not_zero)))
  expect_true(diff(colSums(tab)) == 0)
})
