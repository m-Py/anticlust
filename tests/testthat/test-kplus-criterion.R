
library("anticlust")

context("K-plus anticlustering")

test_that("Test that k-plus criterion is computed correctly", {
  set.seed(123)
  
  init <- initialize_clusters(96, 3, NULL)
  
  df <- schaper2019[, 3:6]
  
  groups1 <- anticlustering(
    df,
    K = init,
    objective = "kplus"
  )
  
  groups2 <- anticlustering(
    cbind(df, squared_from_mean(df)),
    K = init,
    objective = "variance"
  )
  
  expect_true(all(groups1 == groups2))

  # add test with standardize argument
  groups4 <- anticlustering(
    df,
    K = init,
    objective = "kplus",
    standardize = TRUE
  )
  
  groups5 <- anticlustering(
    scale(cbind(df, squared_from_mean(df))),
    K = init,
    objective = "variance",
    standardize = FALSE
  )
  expect_true(all(groups4 == groups5))
  
})
