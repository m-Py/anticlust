
library(anticlust)
library(tinytest)

## Test "vanilla" bin packing

n_batches <- 20
n <- 50
capacities <- rep(10, n_batches)
weights <- sample(5, size = n, replace = TRUE)

opt <- tryCatch(
  anticlust:::optimal_binpacking_(capacities, weights),
  error = function(e) e
)

error <- "simpleError" %in% class(opt)

if (!error) {
  table(opt)
  sums <- tapply(weights, opt, sum)
  sums_all <- rep(0, length(capacities))
  sums_all[as.numeric(names(sums))] <- sums
  # are capacity limits not exceeded?
  expect_true(all(sums_all <= capacities))
}




########################
## Test bin packing as used in must-link anticlustering

must_link_fulfilled <- function(ID, groups_must_link) {
  same <- as.numeric(names(table(ID)[table(ID) > 1]))
  all_good <- c()
  for (i in same) {
    all_good <- c(all_good, all(groups_must_link[ID == i] == groups_must_link[ID == i][1]))
  }
  all(all_good)
}

N <- 120
K <- 10
capacities <- rep(N/K, K)
must_link <- sample(N, replace = TRUE)

start <- Sys.time()
opt <- tryCatch(
  anticlust:::optimal_binpacking_(capacities, table(must_link)),
  error = function(e) e
)
Sys.time() - start

error <- "simpleError" %in% class(opt)

if (!error) {
  table(opt)

  ## This is how we would use optimal bin packing in must-link anticlustering
  # Users cannot explicitly request the optimal method currently, so
  # I emulate how it is done in anticlust
  init <- anticlust:::get_init_assignments(N, ID = must_link, target_groups = capacities, method = "optimal")
  init <- anticlust:::add_unassigned_elements(capacities, init, N, length(capacities))
  # Output of bin packing + filling is as expected?
  expect_true(all(table(init) == capacities))
  # Must-link constraints are met?
  expect_true(must_link_fulfilled(must_link, init))
  
  ## TODO ALWAYS USE TRYCATCH WITH MUST LINK CONSTRAINTS IN CASE THEY CANNOT BE FULFILLED
  # (and then the test cases would fail!)
  groups_must_link <- anticlustering(
    rnorm(N), K = K, must_link = must_link
  )
  expect_true(must_link_fulfilled(must_link, groups_must_link))
  
}
