
library("anticlust")
library(tinytest)
# Function to generate data in test files. Used when numeric imprecision is dangerous
rnd_data_integer <- function(N, M, P = 10, as_factor = FALSE) {
  if (length(P) == 1) {
    data <- data.frame(matrix(sample(P, replace = TRUE, size = N*M), ncol = M))
  } else if (length(P) == M) {
    data <- data.frame(sapply(P, \(x) sample(x, replace = TRUE, size = N)))
  } else {
    stop("Length of P must be 1 or M.")
  }
  if (as_factor) {
    data <- as.data.frame(lapply(data, as.factor))
  }
  colnames(data) <- paste0("X", 1:M)
  data 
}

# all levels of heuristicism work and that exact approach has best objective
conditions <- expand.grid(m = 1:4, p = 2)
for (k in 1:nrow(conditions)) {
  m_features <- conditions[k, "m"]
  p_anticlusters <- conditions[k, "p"]
  n_elements <- p_anticlusters * 5 # n must be multiplier of p
  features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
  ## Traverse through levels of heuristicism
  obj_values <- rep(NA, 4)
  anti_list <- list()
  for (i in c(TRUE, FALSE)) {
    anticlusters <- anticlust:::exact_anticlustering(as.matrix(dist(features)),
                                         p_anticlusters, preclustering = i, NULL)
    anti_list[[i + 1]] <- anticlusters
    # Allow for some numeric imprecision of ILP solver:
    obj_values[i + 1]  <- round(anticlust:::diversity_objective_(anticlusters, features), 10)
  }
  ## Exact solution must have maximum objective
  expect_equal(which.max(obj_values), 1)
}

# Solving ILP works as expected for max dispersion
N <- 30
K <- 3
M <- 5
dat <- rnd_data_integer(N, M)
distances <- dist(dat, method = "manhattan")

ILP <- anticlustering(distances, K = K, objective = "dispersion", method = "ilp")
HEURISTIC <- anticlustering(distances, K = K, objective = "dispersion", method = "local-maximum")

expect_true(dispersion_objective(distances, ILP) >= dispersion_objective(distances, HEURISTIC))
expect_true(all(table(ILP) == N/K))
