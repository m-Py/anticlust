library("anticlust")

context("Greatest common denominator")

# compute input data for clique partitioning
sum_agreements <- function(x) {
  N <- nrow(x)
  M <- ncol(x)
  C_uv <- matrix(ncol = N, nrow = N)
  output_disagreements <- matrix(ncol = N, nrow = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      # Brusco et al. (2009): "In this context, the c_uv values represent the number of attributes on which vertices u and v
      # disagree, minus the number of attributes on which they agree."
      agreements <- sum(x[i, ] == x[j, ], na.rm = TRUE)
      disagreements <- sum(x[i, ] != x[j, ], na.rm = TRUE)
      C_uv[j, i] <- agreements - disagreements
    }
  }
  as.dist(C_uv)
}

test_that("that C implementation of similarity computation is correct", {
  N <- 100
  values <- 5
  M <- 3
  x <- sample(values, size = N * M, replace = TRUE)
  x <- matrix(x, ncol = M)
  results1 <- sum_agreements(x)
  results2 <- compute_cpp_similarities(x)
  expect_true(all(results1 == results2))
})
