#' Number of equal sized partitions
#'
#' @param K How many partitions
#' @param N How many elements
#'
#' @return The number of partitions
#'
#' @export
#' 
#' @examples
#' n_partitions(20, 2)
#'
n_partitions <- function(N, K) {
  if (N %% K != 0) {
    stop("N = ", N, " cannot be evenly split into K = ", K, " sets.")
  }
  x <- rep(N / K, K)
  set_size <- c(sum(x), sum(x) - cumsum(x))
  result   <- 1
  for (i in 1:(length(x)-1)) {
    result <- result * choose(set_size[i], x[i])
  }
  result / factorial(K)
}
