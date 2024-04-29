
library("anticlust")

# argument combinations run through (no preclustering / categories)
methods <- c("exchange", "local-maximum")
conditions <- expand.grid(
  method = methods,
  repetitions = c(1, 5),
  objective = c("variance", "diversity"),
  input = c("distance", "features")
)
illegal_conditions <- conditions$objective == "variance" & conditions$input == "distance"
conditions <- conditions[!illegal_conditions, ]
# for several repetitions, we cannot expect that local maximum always outperforms
# exchange method (latter may have better initial state)
illegal_conditions <- conditions$repetitions > 1 & conditions$method == "local-maximum"
conditions <- conditions[!illegal_conditions, ]

# Set up matrix to store the objective values obtained by different methods
K <- sample(2:6, size = 1)
N <- K * sample(5:10, size = 1)
M <- sample(1:4, size = 1)
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
  conditions[conditions$method != "local-maximum", ], 
  direction = "wide", 
  timevar = "repetitions",
  idvar = colnames(conditions)[c(1, 3, 4)]
)
# test that repetitions improves results
expect_true(all(wide1[, 5] >= wide1[, 4]))

wide2 <- reshape(
  conditions[conditions$repetitions == 1, ], 
  direction = "wide", 
  timevar = "method",
  idvar = colnames(conditions)[2:4]
)
# test that local maximum search tend to improve results
expect_true(all(wide2[, 5] >= wide2[, 4]))


# repeated exchange method works with preclustering / categories
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = 3,
  objective = "variance",
  categories = schaper2019$room,
  repetitions = 2,
  method = "local-maximum"
)
expect_true(all(table(anticlusters, schaper2019$room) == 16))

# Preclustering
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = 4,
  objective = "variance",
  preclustering = TRUE,
  categories = schaper2019$room,
  repetitions = 2,
  method = "local-maximum"
)
expect_true(all(table(anticlusters, schaper2019$room) == 12))

# Test that matched preclusters are assigned to different anticlusters
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = 4,
  objective = "variance",
  preclustering = TRUE,
  method = "local-maximum"
)
matches <- matching(schaper2019[, 3:6], p = 4)
expect_true(all(table(matches, anticlusters) == 1))

anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = 4,
  objective = "variance",
  preclustering = TRUE,
  method = "local-maximum",
  repetitions = 3
)
matches <- matching(schaper2019[, 3:6], p = 4)
expect_true(all(table(matches, anticlusters) == 1))

# Test that matched preclusters are assigned to different anticlusters
# even when a categorical constraint is induced
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  objective = "variance",
  K = 2,
  preclustering = TRUE,
  method = "local-maximum",
  repetitions = 5,
  categories = schaper2019$room
)
matches <- matching(
  schaper2019[, 3:6], 
  p = 2, 
  match_within = schaper2019$room
)
expect_true(all(table(anticlusters, matches)))


# repeated exchange method works with different sized groups
N <- nrow(schaper2019)
n1 <- sample(seq(20, 40, 2), size = 1)
n2 <- sample(seq(20, 40, 2), size = 1)
n3 <- N - n2 - n1
stopifnot((n1 + n2 + n3) == N)

# no repetitions:
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = c(n1, n2, n3),
  categories = schaper2019$room,
  objective = "variance"
)
tab <- table(anticlusters, schaper2019$room)
expect_true(all(tab[, 1] == tab[, 2]))

# only local maximum:
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = c(n1, n2, n3),
  categories = schaper2019$room,
  method = "local-maximum",
  objective = "variance"
)
tab <- table(anticlusters, schaper2019$room)
expect_true(all(tab[, 1] == tab[, 2]))

# only repetitions:
anticlusters <- anticlustering(
  schaper2019[, 3:6],
  K = c(n1, n2, n3),
  categories = schaper2019$room,
  repetitions = 5,
  objective = "variance"
)
tab <- table(anticlusters, schaper2019$room)
expect_true(all(tab[, 1] == tab[, 2]))

# repetitions and local-maximum search:
anticlusters <- anticlustering(
  dist(schaper2019[, 3:6]), # use distances as input at least once here
  K = c(n1, n2, n3),
  categories = schaper2019$room,
  repetitions = 5,
  method = "local-maximum"
)
tab <- table(anticlusters, schaper2019$room)
expect_true(all(tab[, 1] == tab[, 2]))
