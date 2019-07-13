
context("argument preclustering is a vector and not boolean")
library("anticlust")

test_that("Vector input fails for ILP and sampling method", {
  N <- 10
  M <- 2
  features <- matrix(rnorm(N * M), ncol = M)
  expect_error(
    anticlustering(
      features,
      K = 2,
      preclustering = c(rep(1:2, 5)),
      method = "ilp"
    )
  )
  expect_error(
    anticlustering(
      features,
      K = 2,
      preclustering = c(rep(1:2, 5)),
      method = "sampling"
    )
  )
})


test_that("Vector input works for exchange method", {
  N <- 10
  M <- 2
  K <- 2
  features <- matrix(rnorm(N * M), ncol = M)
  for (obj in c("distance", "variance")) {
    preclusters <- get_preclusters(features, NULL, K, TRUE)
    set.seed(123)
    ac1 <- anticlustering(
      features,
      K = K,
      preclustering = preclusters,
      method = "exchange",
      objective = obj
    )
    set.seed(123)
    ac2 <- anticlustering(
      features,
      K = K,
      preclustering = TRUE,
      method = "exchange",
      objective = obj
    )
    ## it is possible that the same assignment has a different return value
    ## because c(1, 1, 2, 2) is the same as c(2, 2, 1, 1). Deal with this:
    if (any(ac1 != ac2)) {
      ## swap 1 and 2 in one of them
      which_one <- ac1 == 1
      ac1[ac1 == 2] <- 1
      ac1[which_one] <- 2
    }
    ## Ensure that the output is the same, regardless of whether preclustering
    ## was computed within or outside of the function
    expect_equal(all(ac1 == ac2), TRUE)
  }
})

