
context("Test repeated exchange methods")
library("anticlust")

test_that("argument combinations run through (no preclustering / categories)", {
  methods <- c("exchange", "local-maximum")
  conditions <- expand.grid(
    method = methods,
    repetitions = c(1, 10),
    objective = c("variance", "diversity"),
    input = c("distance", "features")
  )
  illegal_conditions <- conditions$objective == "variance" & conditions$input == "distance"
  conditions <- conditions[!illegal_conditions, ]
  
  # Set up matrix to store the objective values obtained by different methods
  K <- sample(2:6, size = 1)
  N <- K * sample(7:14, size = 1)
  M <- sample(1:10, size = 1)
  features <- matrix(rnorm(N * M), ncol = M)
  distances <- dist(features)
  objs <- rep(NA, nrow(conditions)) # store objectives
  clusters <- sample(rep_len(1:K, length.out = N))

  for (i in 1:nrow(conditions)) {
    if (conditions$input[i] == "distance") {
      data <- distances
    } else {
      data <- features
    }

    anticlusters <- anticlustering(
      features, 
      K = clusters,
      objective = conditions$objective[i],
      method = conditions$method[i],
      repetitions = conditions$repetitions[i]
    )
    
    if (conditions$objective[i] == "diversity") {
      obj_function <- diversity_objective
    } else {
      obj_function <- variance_objective
    }
    
    objs[i] <- obj_function(data, anticlusters)
  }
  
  # Ensure that repetition exchange method is at least as good as single exchange
  conditions$obj_value <- objs
  # make wide format
  wide1 <- reshape(
    conditions, 
    direction = "wide", 
    timevar = "repetitions",
    idvar = colnames(conditions)[c(1, 3, 4)]
  )
  # test that local maximum search improves results
  expect_true(all(wide1[, 5] >= wide1[, 4]))
  
  wide2 <- reshape(
    conditions, 
    direction = "wide", 
    timevar = "method",
    idvar = colnames(conditions)[2:4]
  )
  # test that repetitions improve results
  expect_true(all(wide2[, 5] >= wide2[, 4]))

})
