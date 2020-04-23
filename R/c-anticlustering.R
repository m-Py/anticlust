
#' C implementation of anticlustering
#' 
#' @param data A N x M data matrix
#' @param clusters An initial assignment of the elements to clusters 
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical constraints. 
#' 
#' @examples
#' # Compare classical anticlustering implementation and C implementation
#' N <- 100
#' M <- 3
#' data <- matrix(rnorm(N * M), ncol = M)
#' K <- 2
#' n_cats <- 2
#' categories <- sample(n_cats, size = N, replace = TRUE)
#' clusters <- anticlust:::categorical_sampling(categories, K)
#' start <- Sys.time()
#' cl1 <- fanticlust(data, clusters, categories)
#' Sys.time() - start
#' table(categories, cl1)
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
fanticlust <- function(data, clusters, categories = NULL) {
  
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
  
  if (argument_exists(categories)) {
    USE_FREQUENCIES <- TRUE
    categories <- merge_into_one_variable(categories) - 1
    N_CATS <- length(unique(categories))
    CAT_frequencies <- table(categories)
  } else {
    USE_FREQUENCIES <- FALSE
    categories <- 0
    N_CATS <- 0
    CAT_frequencies <- 0
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
    as.integer(USE_FREQUENCIES),
    as.integer(N_CATS),
    as.integer(CAT_frequencies),
    as.integer(categories),
    PACKAGE = "anticlust"
  )
  results[["clusters"]] + 1
}