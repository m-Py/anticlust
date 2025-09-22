
library("anticlust")

# Test that k-plus criterion is computed correctly
set.seed(123)

init <- anticlust:::initialize_clusters(96, 3, NULL)

df <- schaper2019[, 3:6]

groups1 <- anticlustering(
  df,
  K = init,
  objective = "kplus"
)

groups2 <- anticlustering(
  kplus_moment_variables(df, T = 2, standardize = FALSE), 
  K = init,
  objective = "variance"
)

expect_true(all(groups1 == groups2))

# add test with standardize argument
groups3 <- anticlustering(
  df,
  K = init,
  objective = "kplus",
  standardize = TRUE
)

groups4 <- anticlustering(
  kplus_moment_variables(df, T = 2), 
  K = init,
  objective = "variance"
)

expect_true(all(groups3 == groups4))
