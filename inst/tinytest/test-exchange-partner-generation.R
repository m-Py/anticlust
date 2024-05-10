

library("anticlust")

# Exchange partners are generated correctly - no restriction
clusters <- c(1, 2, 1, 2, 1, 2)
i <- 1
group_i <- clusters[i]
N <- length(clusters)
categories <- NULL
partners <- anticlust:::get_exchange_partners(clusters, i, categories)
expect_equal(all(partners == c(2, 4, 6)), TRUE)

# Exchange partners are generated correctly - categorical restriction
clusters <- c(1, 2, 1, 2, 1, 2)
i <- 1
group_i <- clusters[i]
N <- length(clusters)
categories <- c(1, 1, 1, 1, 2, 2)
partners <- anticlust:::get_exchange_partners(clusters, i, categories)
expect_equal(all(partners == c(2, 4)), TRUE)

# Exchange partners are generated correctly - preclustering restriction
clusters <- c(1, 2, 1, 2, 1, 2)
i <- 1
group_i <- clusters[i]
N <- length(clusters)
# Add a preclustering restriction (this is now just the same as
# categorical restrictions and does not really make sense here)
preclusters <- c(1, 1, 2, 2, 3, 3)
partners <- anticlust:::get_exchange_partners(clusters, i, preclusters)
expect_equal(partners == 2, TRUE)

# Exchange partners are generated correctly - preclustering and categorical restriction
clusters <- c(1, 2, 1, 2, 1, 2)
i <- 1
group_i <- clusters[i]
N <- length(clusters)
## Use both restrictions
categories <- c(1, 1, 1, 1, 2, 2)
preclusters <- c(1, 2, 2, 1, 3, 3)
constraints <- anticlust:::merge_into_one_variable(cbind(categories, preclusters))
partners <- anticlust:::get_exchange_partners(clusters, i, constraints)
expect_equal(partners == 4, TRUE)

# function that tests that all exchange partners are from the same category
test_idx <- function(idx, categories) {
  if (!is.null(categories)) {
    for (i in 1:length(idx)) {
      category_of_partners <- categories[idx[[i]]]
      # test that all partners have the same category
      expect_true(all(category_of_partners == category_of_partners[1]))
      # test that all partners actually have the category of item `i`
      expect_true(category_of_partners[1] == categories[i])
    }
  }
}

# Exchange partners are generatede correctly for fast k-means method

# generate random categories and data
N <- 1000
M <- 2
k_neighbours <- 10
features <- matrix(rnorm(N * M), ncol = M)
categories <- sample(1:4, size = N, replace = TRUE)

### Using categorical restrictions
# Case 1: restricted number of exchange partners (i.e., nearest neighbour search)

partners <- anticlust:::nearest_neighbours(
  features = features,
  k_neighbours = k_neighbours,
  categories = categories
)
test_idx(partners, categories)
## ensure correct order of the output 
expect_true(all(lengths(partners) == k_neighbours))

# Case 2: no restriction on number of exchange partners
partners <- anticlust:::list_idx_by_category(categories)
test_idx(partners, categories)

### Not using categorical restrictions
# Case 1: restricted number of exchange partners (i.e., nearest neighbour search)
categories <- anticlust:::to_numeric(schaper2019$room)
k_neighbours <- k_neighbours
partners <- anticlust:::nearest_neighbours(
  features = features,
  k_neighbours = k_neighbours,
  categories = NULL
)
expect_true(all(lengths(partners) == k_neighbours))

# remove_self() function works as intended in exchange partner generation
expect_true(
  identical(
    anticlust:::remove_self(list(c(2, 1), c(3, 2), c(4, 3))), list(2, 3, 4)
  )
)
