
library("anticlust")

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
  pattern = "categorical"
)

# no categorical restrictions
expect_error(
  anticlusters <- anticlustering(
    schaper2019[, 3:6],
    K = 4,
    preclustering = TRUE,
    method = "brusco"
  ), 
  pattern = "preclustering"
)

# Test that bicriterion function works as intended -- equal sized groups
bc <- bicriterion_anticlustering(schaper2019[, 3:6], K = 3, R = 20)
tab <- apply(bc, 1, table)
tab <- data.frame(tab)
tab <- tab[order(tab[, 1], decreasing = TRUE), ]
expect_true(all(tab == c(32, 32, 32)))

# Ensure that vector is returned when using "return" argument
expect_true(
  is.vector(bicriterion_anticlustering(schaper2019[, 3:6], K = 3, return = "best-diversity"))
)

expect_true(
  is.vector(bicriterion_anticlustering(schaper2019[, 3:6], K = 3, return = "best-dispersion"))
)

expect_true(is.vector(
  anticlustering(schaper2019[, 3:6], K = 2, method = "brusco")
))

# For dispersion objective
expect_true(is.vector(
  anticlustering(schaper2019[, 3:6], K = 2, method = "brusco", objective = "dispersion")
))

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



## Test return argument

N <- 100
M <- 5
K <- 10
data <- matrix(rnorm(N*M), ncol = M)

set.seed(123)
a <- bicriterion_anticlustering(
  data,
  K = K,
  R = c(10, 10)
)

set.seed(123)
b <- bicriterion_anticlustering(
  data,
  K = K,
  R = c(10, 10)
)

expect_true(identical(a, b))

set.seed(123)
c <- bicriterion_anticlustering(
  data,
  K = K,
  R = c(10, 10),
  return = "best-diversity"
)

set.seed(123)
e <- bicriterion_anticlustering(
  data,
  K = K,
  R = c(10, 10),
  return = "best-dispersion"
)

# Is the 1 solution from the pareto set?
n <- nrow(a)
expect_true(duplicated(rbind(a, c))[n+1])
expect_true(duplicated(rbind(a, e))[n+1])


## Test that computation of the return arguments is correct! 
# Use k-means anticlustering + different sized groups

kmeans_distances <- as.matrix(dist(data)^2)
set.seed(123)
a <- bicriterion_anticlustering(
  kmeans_distances,
  K = c(10, 10, 20, 20, 40),
  R = c(10, 10),
  average_diversity = TRUE
)

set.seed(123)
b <- bicriterion_anticlustering(
  kmeans_distances,
  K = c(10, 10, 20, 20, 40),
  R = c(10, 10),
  average_diversity = TRUE
)

expect_identical(a, b)

set.seed(123)
c <- bicriterion_anticlustering(
  kmeans_distances,
  K = c(10, 10, 20, 20, 40),
  R = c(10, 10),
  average_diversity = TRUE,
  return = "best-average-diversity"
)

n <- nrow(a)
expect_true(duplicated(rbind(a, b))[n+1])

# manually compute objectives:
average_diversities <- c()
diversities <- c()
dispersions <- c()
for (i in 1:nrow(a)) {
  average_diversities[i] <- anticlust:::weighted_diversity_objective_(kmeans_distances, a[i, ], table(a[i, ]))
  diversities[i] <- diversity_objective(kmeans_distances, a[i, ])
}

# Compare with "actual" k-means objective computations
expect_equal(average_diversities, apply(a, 1, variance_objective, x = data)) # numeric imprecision is okay here!

# -> weighted_diversity computation is correct!

# Does BILS function return the partition that maximizes the average diversity?
expect_true(all(a[which.max(average_diversities), ] == c))



# prevent inconsistent calls (average_diversity + return)
expect_error(bicriterion_anticlustering(
  dist(data)^2,
  K = c(10, 10, 20, 20, 40),
  R = c(10, 10),
  average_diversity = FALSE,
  return = "best-average-diversity"
))

expect_error(bicriterion_anticlustering(
  dist(data)^2,
  K = c(10, 10, 20, 20, 40),
  R = c(10, 10),
  average_diversity = TRUE,
  return = "best-diversity"
))
