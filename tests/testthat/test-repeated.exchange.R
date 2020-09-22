
context("Test repeated exchange methods")
library("anticlust")

test_that("all argument combinations run through", {
  methods <- c("exchange", "local-maximum")
  conditions <- expand.grid(
    preclustering = c(TRUE, FALSE),
    method = methods,
    repetitions = c(1, 10),
    objective = c("variance", "diversity"),
    input = c("distance", "features")
  )
  illegal_conditions <- conditions$objective == "variance" & conditions$input == "distance"
  conditions <- conditions[!illegal_conditions, ]
  
  # Set up matrix to store the objective values obtained by different methods
  K <- sample(2:6, size = 1)
  N <- K * sample(10:20, size = 1)
  M <- sample(1:10, size = 1)
  features <- matrix(rnorm(N * M), ncol = M)
  distances <- dist(features)
  objs <- rep(NA, nrow(conditions)) # store objectives
  initial_clusters <- sample(rep_len(1:K, length.out = N))
  preclusters <- matching(features, K)
  initial_clusters_preclustered <- initialize_clusters(N, K, preclusters)
  
  for (i in 1:nrow(conditions)) {
    if (conditions$input[i] == "distance") {
      data <- distances
    } else {
      data <- features
    }
    
    if (conditions$preclustering[i] == TRUE) {
      categories <- preclusters
      clusters <- initial_clusters_preclustered
    } else {
      categories <- NULL
      clusters <- initial_clusters
    }
    
    anticlusters <- anticlustering(
      features, 
      K = clusters,
      objective = conditions$objective[i],
      method = conditions$method[i],
      categories = preclusters,
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
    idvar = colnames(conditions)[c(1, 2, 4, 5)]
  )
  wide1
  wide2 <- reshape(
    conditions, 
    direction = "wide", 
    timevar = "method",
    idvar = colnames(conditions)[c(1, 3, 4, 5)]
  )
  wide2
  
  
})
