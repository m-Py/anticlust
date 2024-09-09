
n_batches <- 5
n <- 6
capacities <- rep(5, n_batches)
weights <- sample(5, size = n, replace = TRUE)

opt <- tryCatch(
  optimal_binpacking_(capacities, weights),
  error = function(e) e
)

error <- "simpleError" %in% class(opt)

if (!error) {
  sums <- tapply(weights, opt, sum)
  sums_all <- rep(0, length(capacities))
  sums_all[as.numeric(names(sums))] <- sums
  expect_true(all(sums_all <= capacities))
}
