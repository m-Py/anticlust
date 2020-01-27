
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The `testthat` package is not available, therefore `anticlust` cannot be tested.")
}

library("anticlust")
library("testthat")

test_check("anticlust")

