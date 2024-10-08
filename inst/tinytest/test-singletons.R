
library(tinytest)
library(anticlust)

# Since version (0.8.7) it is possible that a singleton cluster exists in anticlustering().
# This test files tries to make sure that this does not cause any problems (prior to 
# that, requesting a singleton cluster just threw an error because it can also 
# be left out completely from anticlustering. However, with the implementation of 
# must-link constraints, it is now possible that we "accidentally" have a singleton 
# without intention of the user cluster)

features <- schaper2019[, 3:6]
cluster_sizes <- c(1, 1, 1, 49, nrow(features) - (49+1+1+1))
tt <- anticlustering(
  features,
  K = cluster_sizes
)
expect_true(all(sort(table(tt)) == sort(cluster_sizes)))

tt <- anticlustering(
  features,
  K = cluster_sizes,
  objective = "dispersion"
)
expect_true(all(sort(table(tt)) == sort(cluster_sizes)))

tt <- anticlustering(
  features,
  K = cluster_sizes,
  objective = "variance"
)
expect_true(all(sort(table(tt)) == sort(cluster_sizes)))

tt <- anticlustering(
  features,
  K = cluster_sizes,
  objective = "kplus"
)
expect_true(all(sort(table(tt)) == sort(cluster_sizes)))



# Ensure that diversity_objective is computed correctly when there are singleton clusters

expect_true(
  sum(anticlust:::diversity_objective_by_group(tt, features)) == diversity_objective(features, tt)
)

#############

# C implementation has same result as R implementation when there are singletons?
init <- sample(tt)

tt <- anticlustering(
  features,
  K = init
)

uu <- anticlustering(
  features, 
  K = init,
  objective = diversity_objective
)

expect_true(all(tt == uu))


### Do these tests for other objectives:

# Variance

init <- sample(tt)

tt <- anticlustering(
  features,
  K = init,
  objective = "variance"
)

uu <- anticlustering(
  features, 
  K = init,
  objective = variance_objective
)

expect_true(all(tt == uu))

# Dispersion

init <- sample(tt)

tt <- anticlustering(
  features,
  K = init,
  objective = "dispersion"
)

uu <- anticlustering(
  features, 
  K = init,
  objective = dispersion_objective
)

expect_true(all(tt == uu))


# - Test for singleton created via must-link constraint
N <- nrow(features)
K <- 12
must_link <- rep(NA, N)
must_link[1:(N/K)] <- 1
cl <- anticlustering(features, K = K, must_link = must_link)
expect_true(all(cl[1] == cl[1:(N/K)]))
expect_true(all(cl[1] != cl[((N/K)+1):N]))

# now make it impossible to fulfill the constraints:
must_link[1:((N/K)+1)] <- 1
expect_error(
  anticlustering(features, K = K, must_link = must_link),
  pattern = "must-link"
)
