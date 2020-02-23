

context("Generate exchange partners")
library("anticlust")

test_that("Exchange partners are generated correctly - no restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- NULL
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4, 6)), TRUE)
})

test_that("Exchange partners are generated correctly - categorical restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  categories <- c(1, 1, 1, 1, 2, 2)
  partners <- get_exchange_partners(clusters, i, categories)
  expect_equal(all(partners == c(2, 4)), TRUE)
})

test_that("Exchange partners are generated correctly - preclustering restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  # Add a preclustering restriction (this is now just the same as
  # categorical restrictions and does not really make sense here)
  preclusters <- c(1, 1, 2, 2, 3, 3)
  partners <- get_exchange_partners(clusters, i, preclusters)
  expect_equal(partners == 2, TRUE)
})

test_that("Exchange partners are generated correctly - preclustering and categorical restriction", {
  clusters <- c(1, 2, 1, 2, 1, 2)
  i <- 1
  group_i <- clusters[i]
  N <- length(clusters)
  ## Use both restrictions
  categories <- c(1, 1, 1, 1, 2, 2)
  preclusters <- c(1, 2, 2, 1, 3, 3)
  constraints <- merge_into_one_variable(cbind(categories, preclusters))
  partners <- get_exchange_partners(clusters, i, constraints)
  expect_equal(partners == 4, TRUE)
})

# function that tests that all exchange partners are from the same category
test_idx <- function(idx, categories) {
  if (argument_exists(categories)) {
    for (i in 1:length(idx)) {
      category_of_partners <- categories[idx[[i]]]
      # test that all partners have the same category
      expect_true(all(category_of_partners == category_of_partners[1]))
      # test that all partners actually have the category of item `i`
      expect_true(category_of_partners[1] == categories[i])
    }
  }
}

test_that("Exchange partners are generatede correctly for fast k-means method", {

  # generate random categories and data
  N <- 1000
  M <- 2
  k_neighbours <- 10
  features <- matrix(rnorm(N * M), ncol = M)
  categories <- sample(1:4, size = N, replace = TRUE)
  
  ### Using categorical restrictions
  # Case 1: restricted number of exchange partners (i.e., nearest neighbour search)
  
  partners <- all_exchange_partners(
    features = features,
    k_neighbours = k_neighbours,
    categories = categories
  )
  test_idx(partners, categories)
  ## ensure correct order of the output 
  expect_true(all(sapply(partners, function(x) x[1]) == 1:N))
  
  # Case 2: no restriction on number of exchange partners
  partners <- all_exchange_partners(
    features = features,
    k_neighbours = Inf,
    categories = categories
  )
  test_idx(partners, categories)
  
  ### Not using categorical restrictions
  # Case 1: restricted number of exchange partners (i.e., nearest neighbour search)
  categories <- to_numeric(schaper2019$room)
  k_neighbours <- k_neighbours
  partners <- all_exchange_partners(
    features = features,
    k_neighbours = k_neighbours,
    categories = NULL
  )
  expect_true(all(sapply(partners, length) == k_neighbours + 1))
  ## ensure correct order of the output 
  expect_true(all(sapply(partners, function(x) x[1]) == 1:N))
  
  # Case 2: no restriction on number of exchange partners
  partners <- all_exchange_partners(
    features = features,
    k_neighbours = Inf,
    categories = NULL
  )
  expect_true(all(sapply(partners, length) == N))
})


# test for external data set. This test used to fail due to 
# categories very very few members. useful test on actual data
test_that("Exchange partners are generatede correctly for fast k-means method", {
  skip_on_cran()
  temp <- tempfile()
  download.file("http://openpsychometrics.org/_rawdata/NPI.zip", temp)
  npi <- read.csv(unz(temp, "NPI/data.csv")) # reads the data into R
  unlink(temp)
  rm(temp)
  
  temp <- npi[, paste0("Q", 1:40)]
  temp[temp == 0] <- NA # 0 = missing value in this data set
  npi <- npi[complete.cases(temp), ]
  rm(temp)
  
  # The following vector encodes the key (i.e., narcissistic response)
  # by item:
  keys <- c(1, 1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2,
            1, 2, 2, 1, 1, 2, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 2)
  
  # Use for-loop to score all 40 items using the key
  for (i in 1:40) {
    colname <- paste0("Q", i)
    npi[[paste0("score_", colname)]] <- ifelse(npi[[colname]] == keys[i], 1, 0)
  }
  
  item_responses <- subset(npi, select = score_Q1:score_Q40)
  
  partners <- all_exchange_partners(
    features = item_responses,
    k_neighbours = 5,
    categories = to_numeric(npi$gender)
  )
  
  test_idx(partners, to_numeric(npi$gender))
})
