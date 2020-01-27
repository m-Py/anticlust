
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The `testthat` package is not available, therefore `anticlust` cannot be tested.")
}

library("anticlust")

test_check("anticlust")

