
library(tinytest)
library(anticlust)

# Since version (0.8.7) it is possible that a singleton cluster exists in anticlustering
# make sure, this does not cause any problems (prior to that, requesting a singleton
# cluster just threw an error because it can also be left out completely from anticlustering.
# but with the implementation of must-link constraints, it is possible that we "accidentally"
# have a singleton cluster)
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


### TODO: Do these tests for other objectives:

# - Ensure that objective is computed correctly when there are singleton clusters
# - C implementation has same result as R implementation when there are singletons
