
library("anticlust")

conditions <- expand.grid(m = 1:4, p = 2)
for (k in 1:nrow(conditions)) {
  m_features <- conditions[k, "m"]
  K <- conditions[k, "p"]
  n_elements <- K * 5 # n must be multiplier of p
  features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
  distances <- dist(features)
  
  ac_feat <- anticlustering(features, K = K, preclustering = FALSE,
                            method = "ilp")
  ac_dist <- anticlustering(distances, K = K,
                            preclustering = FALSE,
                            method = "ilp")
  expect_equal(diversity_objective(distances, ac_feat),
               diversity_objective(features, ac_feat))
  expect_equal(diversity_objective(distances, ac_dist),
               diversity_objective(features, ac_dist))
  expect_equal(diversity_objective(distances, ac_feat),
               diversity_objective(distances, ac_dist))
}


conditions <- expand.grid(m = 1:4, p = 2)
for (k in 1:nrow(conditions)) {
  m_features <- conditions[k, "m"]
  K <- conditions[k, "p"]
  n_elements <- K * 5 # n must be multiplier of p
  features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
  distances <- dist(features)
  
  ac_feat <- anticlustering(features, K = K, preclustering = TRUE,
                            method = "ilp")
  ac_dist <- anticlustering(distances, K = K,
                            preclustering = TRUE,
                            method = "ilp")
  expect_equal(diversity_objective(distances, ac_feat),
               diversity_objective(features, ac_feat))
  expect_equal(diversity_objective(distances, ac_dist),
               diversity_objective(features, ac_dist))
  expect_equal(diversity_objective(distances, ac_feat),
               diversity_objective(distances, ac_dist))
  
  # ensure that preclusters are balanced between anticlusters
  preclusters <- balanced_clustering(features, K = n_elements / K, method = "ilp")
  expect_true(all(table(ac_feat, preclusters) == 1))
}


conditions <- expand.grid(m = 1:4, p = 2:4)
for (k in 1:nrow(conditions)) {
  m_features <- conditions[k, "m"]
  K <- conditions[k, "p"]
  n_elements <- K * 5 # n must be multiplier of p
  features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
  distances <- dist(features)
  
  ## Use a fixed seed to compare the random sampling method based on
  ## features and distance input
  rnd_seed <- sample(10000, size = 1)
  
  set.seed(rnd_seed)
  ac_feat <- anticlustering(features, K = K)
  set.seed(rnd_seed)
  ac_dist <- anticlustering(distances, K = K)
  
  expect_equal(diversity_objective(distances, ac_feat),
               diversity_objective(features, ac_feat))
  expect_equal(diversity_objective(distances, ac_dist),
               diversity_objective(features, ac_dist))
  expect_equal(diversity_objective(distances, ac_feat),
               diversity_objective(distances, ac_dist))
}

# At the end, just test that the input works without error when preclustering is enabled
anticlustering(features, K = K, preclustering = TRUE)
anticlustering(distances, K = K, preclustering = TRUE)
