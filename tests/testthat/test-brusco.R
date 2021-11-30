
context("Brusco et al. algorithm")
library("anticlust")

test_that("input and output work expectedly for Brusco algorithm", {
  
  # Test standard usage: K = number of groups
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = 3,
    method = "brusco"
  )
  expect_true(all(table(anticlusters) == c(32, 32, 32)))
  
  # Different group sizes 
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = c(48, 24, 24),
    method = "brusco"
  )
  expect_true(all(table(anticlusters) == c(48, 24, 24)))

  # no categorical restrictions
  expect_error(
    anticlusters <- anticlustering(
      schaper2019[, 3:6],
      K = 4,
      categories = schaper2019$room,
      method = "brusco"
    ), 
    regexp = "categorical"
  )
  
  # no categorical restrictions
  expect_error(
    anticlusters <- anticlustering(
      schaper2019[, 3:6],
      K = 4,
      preclustering = TRUE,
      method = "brusco"
    ), 
    regexp = "preclustering"
  )
  
  # Test that bicriterion function works as intended -- equal sized groups
  bc <- bicriterion_anticlustering(schaper2019[, 3:6], K = 3, R = 20)
  tab <- apply(bc, 1, table)
  tab <- data.frame(tab)
  tab <- tab[order(tab$X1, decreasing = TRUE), ]
  expect_true(all(tab == c(32, 32, 32)))
  
  
  # Test that bicriterion function works as intended -- unequal sized groups
  bc <- bicriterion_anticlustering(schaper2019[, 3:6], K = c(12, 12, 12, 60), R = 20)
  tab <- apply(bc, 1, table)
  tab <- data.frame(tab)
  tab <- tab[order(tab$X1, decreasing = TRUE), ]
  expect_true(all(tab == c(60, 12, 12, 12)))
  
  # Ensure that random seeds work with bicriterion algorithm (this is important
  # because the algorithm uses random number generation in C, should be 
  # reproducible from R!)
  set.seed(1)
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    method = "brusco"
  )
  set.seed(1)
  anticlusters2 <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    method = "brusco"
  )
  expect_true(all(anticlusters == anticlusters2))
  
  # Other seed = different results
  set.seed(3)
  anticlusters3 <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    method = "brusco"
  )
  expect_true(!all(anticlusters == anticlusters3))
  
  # Same with the dispersion criterion
  # Ensure that random seeds work with bicriterion algorithm
  set.seed(1)
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    method = "brusco",
    objective = "dispersion"
  )
  set.seed(1)
  anticlusters2 <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    method = "brusco",
    objective = "dispersion"
  )
  expect_true(all(anticlusters == anticlusters2))
  
  # Other seed = different results
  set.seed(3)
  anticlusters3 <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    method = "brusco",
    objective = "dispersion"
  )
  expect_true(!all(anticlusters == anticlusters3))
  
})
