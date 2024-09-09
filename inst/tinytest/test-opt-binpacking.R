
library(anticlust)
library(tinytest)

## Test "vanilla" bin packing

n_batches <- 5
n <- 6
capacities <- rep(5, n_batches)
weights <- sample(5, size = n, replace = TRUE)

opt <- tryCatch(
  anticlust:::optimal_binpacking_(capacities, weights),
  error = function(e) e
)
table(opt)

error <- "simpleError" %in% class(opt)

if (!error) {
  sums <- tapply(weights, opt, sum)
  sums_all <- rep(0, length(capacities))
  sums_all[as.numeric(names(sums))] <- sums
  expect_true(all(sums_all <= capacities))
}

## Test bin packing as used in must-link anticlustering

N <- 120
K <- 10
capacities <- rep(N/K, K)
must_link <- sample(N, replace = TRUE)

start <- Sys.time()
opt <- tryCatch(
  anticlust:::optimal_binpacking_(capacities, table(must_link), solver = "lpSolve"),
  error = function(e) e
)
Sys.time() - start
table(opt)

## TODO ALWAYS USE TRYCATCH WITH MUST LINK CONSTRAINTS IN CASE THEY CANNOT BE FULFILLED
# (and then the test cases would fail!)
groups_must_link <- anticlustering(
  rnorm(N), K = K, must_link = must_link
)

ID <- must_link
same <- as.numeric(names(table(ID)[table(ID) > 1]))
for (i in same) {
  expect_true(all(groups_must_link[ID == i] == groups_must_link[ID == i][1]))
}
