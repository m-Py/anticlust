
# Function to generate data in test files. Do not use normal data due to numeric imprecision

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

