#' Get the next permutation 
#'
#' Each permutation is computed on basis of a passed permutation where
#' lexicographic ordering of the permutations is used to determine the
#' next "higher" permutation. This is an adaption of the
#' `next_permutation` function in C++.
#'
#' @param permutation A vector of elements.
#'
#' @return The next higher permutation of the elements in vector
#'     `permutation` with regard to its lexicographic ordering.
#'
#' 
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @noRd
#' 
next_permutation <- function(permutation) {
  n    <- length(permutation)
  last <- permutation[n]
  i    <- n
  while(last <= permutation[i-1]) {
    last <- permutation[i-1]
    i    <- i - 1
    ## if lexicographic order is already at the maximum:
    if (i-1 == 0) return(sort(permutation))
  }
  ## this algorithm divides the input in a head and a tail; the tail
  ## is monotonically decreasing
  head <- permutation[1:(i-1)]
  tail <- permutation[i:length(permutation)]
  ## which element in the tail is the smallest element that is larger
  ## than the last element in the head?
  larger_values  <- tail[tail > head[length(head)]]
  ## last element of the head:
  final_head <- head[length(head)]
  ## replace last element of head by smallest larger value in tail
  head[length(head)] <- min(larger_values)  
  ## replace smallest larger value in tail by final head element
  tail[max(which(tail == min(larger_values)))] <- final_head
  ## reverse tail before returning
  return(c(head, rev(tail)))
}
