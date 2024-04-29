
context("Categorical constraints")
library("anticlust")

test_that("categorical constraints are met for one categorie (vector input)", {
  N <- 60
  M <- 2
  features <- matrix(rnorm(N * M), ncol = M)
  ## iterate over number of categories & number of anticlusters
  for (K in 2:4) {
    for (C in 2:4) {
      categories <- sample(rep_len(1:C, N))
      ac <- categorical_sampling(categories = categories, K = K)
      tab <- table(categories, ac)
      ## At most 1 deviation between categories in anticlusters
      expect_equal(all(abs(tab - tab[[1]]) <= 1), TRUE)
    }
  }
})


## The following test assumes that the categorical variables are
## balanced (otherwise testing is hard)

test_that("categorical constraints are met for two categories (data.frame or matrix input)", {
  M <- 2
  ## iterate over number of categories & number of anticlusters
  for (K in 2:3) {
    for (C1 in 2:3) {
      for (C2 in 2:3) {
        ## 1. Choose appropriate N that allows for a balanced assignment
        N <- (K * C1 * C2)^2
        features <- matrix(rnorm(N * M), ncol = M)
        ## 2. Ensure that the categories are actually balanced
        categories1 <- sort(rep_len(1:C1, N))
        categories2 <- categorical_sampling(categories1, K)
        frame_together <- ifelse(sample(2, size = 1) <= 2, data.frame, cbind)
        categories <- frame_together(categories1, categories2)
        ## Random order to cath potential problem in the implementation
        ## that might be due to a reliance on sorted input (We don't
        ## want that)
        categories <- categories[sample(nrow(categories)), ]
        vectorized_categories <- factor(do.call(paste0, as.list(categories)))
        ## Critical test: did it work to create balanced groups?
        tab <- table(vectorized_categories)
        expect_equal(all(tab == tab[[1]]), TRUE)


        ## 3. Are the categories balanced across anticlusters?
        ac <- categorical_sampling(K = K, categories = merge_into_one_variable(categories))
        tab1 <- table(categories[, 1], ac)
        tab2 <- table(categories[, 2], ac)
        tab3 <- table(categories[, 1], categories[, 2], ac)
        ## At most 1 deviation between categories in anticlusters
        expect_equal(all(abs(tab1 - tab1[[1]]) <= 1), TRUE)
        expect_equal(all(abs(tab2 - tab2[[1]]) <= 1), TRUE)
        expect_equal(all(abs(tab3 - tab3[[1]]) <= 1), TRUE)
      }
    }
  }
})


test_that("argument categories can be used to enforce preclustering constraints correctly", {
  N <- 60
  M <- 2
  features <- matrix(rnorm(N * M), ncol = M)
  ## iterate over number of categories & number of anticlusters
  for (K in 2:4) {
    preclusters <- balanced_clustering(features, K = N / K)
    ac <- categorical_sampling(K = K, categories = preclusters)
    tab <- table(ac, preclusters)
    expect_equal(all(tab == 1), TRUE)
  }
})

test_that("balanced categorical output for all anticlustering combinations", {
  features <- schaper2019[, 3:6]
  K <- 3
  categories <- schaper2019$room
  # Anticluster editing
  ac <- anticlustering(
    features,
    K = K,
    categories = categories
  )
  expect_true(all(table(ac, categories) == 16))
  
  # Anticluster editing, preclustering
  ac <- anticlustering(
    features,
    K = K,
    categories = categories,
    preclustering = TRUE
  )
  expect_true(all(table(ac, categories) == 16))
  
  # Anticluster editing, distance input
  ac <- anticlustering(
    dist(features),
    K = K,
    categories = categories
  )
  expect_true(all(table(ac, categories) == 16))
  
  # Anticluster editing, distance input, preclustering
  ac <- anticlustering(
    dist(features),
    K = K,
    categories = categories,
    preclustering = TRUE
  )
  expect_true(all(table(ac, categories) == 16))
  
  # K-means anticlustering 
  ac <- anticlustering(
    features,
    K = K,
    categories = categories
  )
  expect_true(all(table(ac, categories) == 16))
  
  # K-means anticlustering, preclustering
  ac <- anticlustering(
    features,
    K = K,
    categories = categories,
    preclustering = TRUE
  )
  expect_true(all(table(ac, categories) == 16))
  
  # Fast k-means anticlustering 
  ac <- fast_anticlustering(
    features,
    K = K,
    categories = categories
  )
  expect_true(all(table(ac, categories) == 16))
  
  # Fast k-means anticlustering, reduced number of exchange partners
  ac <- fast_anticlustering(
    features,
    K = K,
    categories = categories,
    k_neighbours = 10
  )
  expect_true(all(table(ac, categories) == 16))
})

test_that("Categories argument works with kplus_anticlustering()", {
  # kplus_anticlustering must work as well 
  features <- schaper2019[, 3:6]
  anticlusters <- kplus_anticlustering(features, K = 3, categories = schaper2019$room)

  expect_true(all(table(schaper2019$room, anticlusters) == 16))
  
  # does preclustering work at the same time as well?
  
  anticlusters <- kplus_anticlustering(features, K = 3, categories = schaper2019$room, preclustering = TRUE)
  matches <- matching(features, p = 3, match_within = schaper2019$room)
  expect_true(all(table(schaper2019$room, anticlusters) == 16))
  expect_true(all(table(matches, anticlusters) == 1))

})
