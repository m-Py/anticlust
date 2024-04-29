
context("Standardize argument")
library("anticlust")

test_that("standardize argument only works for intended input", {
  expect_error(
    anticlustering(
      rnorm(10),
      K = 2,
      standardize = 123
    )
  )
  
  expect_error(
    anticlustering(
      rnorm(10),
      K = 2,
      standardize = "foof"
    )
  )
  
  expect_error(
    anticlustering(
      rnorm(10),
      K = 2,
      standardize = c(TRUE, FALSE)
    )
  )
  
  expect_error(
    anticlustering(
      rnorm(10),
      K = 2,
      standardize = NA
    )
  )
  
  expect_true(length(
    anticlustering(
      rnorm(10),
      K = 2,
      standardize = TRUE
    )) == 10)
  
  expect_true(length(
    anticlustering(
      rnorm(10),
      K = 2,
      standardize = FALSE
    )) == 10)

})
