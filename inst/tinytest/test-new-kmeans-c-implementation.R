library("anticlust")
library("RANN")

# New k-means C implementation yields same results as other implementations
  
  set.seed(123)
  M <- 5
  N <- 180
  K <- 4
  features <- matrix(rnorm(N*M), ncol = M)
  init <- sample(rep_len(1:K, nrow(features)))
  

  # fast_anticlustering() does the same as anticlustering()?
  ac_exchangeC <- fast_anticlustering(
    features, K = init, exchange_partners = rep(list(1:N), N) # -> everyone is exchange partner
  )
  ac_exchangeC2 <- anticlustering(features, K = init, objective = "variance")

  expect_true(all(ac_exchangeC == ac_exchangeC2))
  
  # new C fast_anticlustering() implementation does the same as old R implementation?
  ac_exchangeC <- fast_anticlustering(features, K = init, k_neighbours = 20)
  ac_exchangeR <- anticlust:::fast_exchange_(features, init, anticlust:::nearest_neighbours(features, 20, NULL))
  expect_true(all(ac_exchangeC == ac_exchangeR))
  
  # Also test with categorical restrictions
  features <- as.matrix(schaper2019[, 3:6])
  categories <- schaper2019$room
  N <- nrow(features)
  ac_exchangeC <- fast_anticlustering(features, K = 4, k_neighbours = N-1, categories = categories)
  ac_exchangeC2 <- anticlustering(features, K = 4, objective = "variance", categories = categories)
  expect_true(all(table(ac_exchangeC, categories)  == table(ac_exchangeC2, categories)))
  
  # Use reduced exchange partners
  init <- anticlust:::initialize_clusters(N, 3, categories)
  ac_exchangeC <- fast_anticlustering(features, K = init, k_neighbours = 20, categories = categories)
  ac_exchangeR <- anticlust:::fast_exchange_(features, init, anticlust:::nearest_neighbours(features, 20, anticlust:::to_numeric(categories)))
  expect_true(all(ac_exchangeC == ac_exchangeR))
  
  # What if `k_neighbours` is greater than the number of elements in the group with fewest members?
  categories <- anticlust:::merge_into_one_variable(cbind(schaper2019$syllables, schaper2019$room))
  init <- categorical_sampling(categories, K = 3)
  ac_exchangeC <- fast_anticlustering(features, K = init, k_neighbours = 10, categories = categories)
  ac_exchangeR <- anticlust:::fast_exchange_(features, init, anticlust:::nearest_neighbours(features, 10, categories))
  expect_true(all(table(ac_exchangeC, categories) == table(ac_exchangeR, categories)))
  
  # test custom exchange partners
  N <- nrow(features)
  # random exchange partners
  k_neighbours <- 10
  exchange_partners <- lapply(rep(N, N), function(x) sample(1:x)[1:k_neighbours])
  # R and C implementation the same?
  ac_exchangeC <- fast_anticlustering(features, K = init, exchange_partners = exchange_partners)
  ac_exchangeR <- anticlust:::fast_exchange_(features, init, exchange_partners)
  expect_true(all(ac_exchangeC == ac_exchangeR))
  
  # Remove some exchange partners to test if an unequal number of exchange partners works
  exchange_partners[[1]] <- exchange_partners[[1]][-(1:3)]
  exchange_partners[[10]] <- exchange_partners[[10]][-2]
  exchange_partners[[23]] <- exchange_partners[[23]][-(10:7)]
  ac_exchangeC <- fast_anticlustering(features, K = init, exchange_partners = exchange_partners)
  ac_exchangeR <- anticlust:::fast_exchange_(features, init, exchange_partners)
  expect_true(all(ac_exchangeC == ac_exchangeR))
  
  # Does the argument `exchange_partners` really do the same as `k_neighbours`?
  nn_exchange_partners <- RANN::nn2(features, k = k_neighbours + 1)$nn.idx #+1 because the element itself is also returned by nn2
  ac1 <- fast_anticlustering(features, K = init, exchange_partners = anticlust:::matrix_to_list(nn_exchange_partners))
  ac2 <- fast_anticlustering(features, K = init, k_neighbours = k_neighbours)
  expect_true(all(ac1 == ac2))
  
  # ensure that using exchange_partners does the same as not using it, when combining with new function generate_exchange_partners()
  ac1 <- fast_anticlustering(features, K = init, exchange_partners = generate_exchange_partners(5, features = features, method = "RANN"))
  ac2 <- fast_anticlustering(features, K = init, k_neighbours = 5)
  expect_true(all(ac1 == ac2))
  
  # again, using a different data set
  features <- iris[,-5]
  init <- anticlust:::initialize_clusters(K = 75, N = nrow(features), NULL)
  n_exchange_partners <- 10
  groups1 <- fast_anticlustering(
    features,
    K = init,
    exchange_partners = generate_exchange_partners(n_exchange_partners, features = features, method = "RANN")
  )
  # compare with using nearest neighbours as exchange partners:
  groups2 <- fast_anticlustering(features, K = init, k_neighbours = n_exchange_partners)
  expect_true(all(groups1 == groups2))
  
  ## Test that new default for k_neighbours works (and if N < 20, it still works!)
  expect_length(fast_anticlustering(1:6, K = 2), 6)
  expect_length(fast_anticlustering(1:20, K = 2), 20)
  expect_length(fast_anticlustering(1:21, K = 2), 21)
