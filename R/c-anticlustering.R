
#' C implementation of anticlustering
#' 
#' @param data A N x M data matrix
#' @param clusters An initial assignment of the elements to clusters 
#' 
#' @examples
#' # Compare classical anticlustering implementation and C implementation
#' N <- 500
#' M <- 3
#' data <- matrix(rnorm(N * M), ncol = M)
#' K <- 3
#' clusters <- sample(rep_len(1:K, N))
#' start <- Sys.time()
#' cl1 <- fanticlust(data, clusters)
#' Sys.time() - start
#'
#' start <- Sys.time()
#' cl2 <- anticlustering(
#'   data,
#'   K = clusters,
#'   objective = "variance"
#' )
#' Sys.time() - start
#' 
#' all(cl1 == cl2) # Same algorithm in C and R
#' 
#' 
#' @useDynLib anticlust c_anticlustering
#' 
#' @export
#' 
fanticlust <- function(data, clusters) {
  
  data <- as.matrix(data)
  if (mode(data) != "numeric") {
    stop("The data must be numeric.")
  }
  
  # Ensure that `clusters` is a vector of integers between 0 and K-1
  clusters <- to_numeric(clusters) - 1
  # As convenience arguments for C, also pass 
  # a. the number of clusters (`K`)
  # b. how often each cluster occurs (`frequencies`)
  # c. the number of rows (`N`) and columns (`M`) in `data`
  K <- length(unique(clusters))
  frequencies <- table(clusters)
  N <- NROW(data)
  M <- NCOL(data)
  
  # Are the arguments compatible?
  if (N != length(clusters)) {
    stop("The two arguments imply different length of the data.")
  }
  
  # Each cluster must occur more than once
  if (any(frequencies <= 1)) {
    stop("No clusters with only one member allowed.")
  }
  
  # Call C implementation of anticlustering
  results <- .C(
    "c_anticlustering", 
    as.double(data),
    as.integer(N),
    as.integer(M),
    as.integer(K),
    as.integer(frequencies),
    clusters = as.integer(clusters),
    PACKAGE = "anticlust"
  )
  results[["clusters"]] + 1
}
