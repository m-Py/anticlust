  #' Solve anticlsuter editing by using the BILS algorithm by Brusco implemented in c
#' 
#' @param data A N x N dissimilarity matrix or N x M features matrix.
#' @param G The number of clusters 
#' @param R The desired number of restarts for the algorithm
#' @param (optional) W list of possible weights for dispersion in the bicriterion (0 <= W <= 1)
#' @param (optional) Xi a vector for the neighborhood size range [xi1,xi2] in percent 
#' (for example c(0.05, 0.1) can be interpreted as 5%-10%)
#' 
#' @return list of useful partitions (paretoset)
#' 
#' @useDynLib anticlust
#' 
#' @export
  


bicriterion_iterated_local_search_call <- function(data, G, R, W = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 0.99, 0.999, 0.999999), Xi = c(0.05,0.1)) {
  
  distances <- convert_to_distances(data) 
  N <- NROW(distances)
  WL <- length(W)
  
  checkweights(W);
  checkneighborhood(Xi);
  
  if((N %% G) != 0){
    stop("The number of elements must be divisble by the number of clusters.")
  }
  
  upper_bound = 500 # limits number of resulting partitions 
  result_matrix = matrix(data = -1, nrow = upper_bound , ncol = N) # create empty matrix for results to use in c
  
  # Call C function
  bicriterion <- .C(
    "bicriterion_iterated_local_search_call",
    as.double(distances),
    as.integer(N),
    as.integer(G),
    as.integer(R),
    as.integer(upper_bound),
    as.integer(WL),
    as.double(W),
    as.double(Xi),
    result = double(length(result_matrix)),
    PACKAGE = "anticlust" # important to call C
  )

  #c returns the list of partitions as one vector, that we turn back into a matrix
  result_matrix <- vector_to_matrix(bicriterion[["result"]],result_matrix, upper_bound, N)

  return(result_matrix)
}


checkweights <- function(weights){
  for(i in weights){
    if(i < 0 | i > 1){
      stop("All weights must be a percentage. Only values between 0 and 1 are allowed.")
    }
  }  
}

checkneighborhood <- function(Xi){
  for(i in Xi){
    if(i < 0 | i > 1){
      stop("Both neighorhoodindexes must be a percentage. Only values between 0 and 1 are allowed.")
    }
  } 
  if(Xi[1] > Xi[2]){
    stop("First neighborhoodpercentage needs to be smaller than the second.")
  }
}


#fill matrix with results and delete unnecessary columns 
vector_to_matrix <- function(vector, matrix, upper_bound, size){
  row <- 1
  column <- 1
  pos <- 1
  while(vector[pos] != -1 & row <= upper_bound ){
    while(column<= size){
      matrix[row, column] <- vector[pos]
      pos <- pos + 1
      column <- column + 1
    }
    column <- 1
    row <- row + 1
  }
  
  return(matrix[1:row-1,1:size])
}


plot_partition_matrix <- function(matrix, data, restarts){
  
  distances <- convert_to_distances(data)
  
  if(NCOL(matrix) == 1){
    diversity <- get_diversity(matrix, distances)
    dispersion <- get_dispersion(matrix, distances)
    matrix <- cbind(diversity, dispersion, restarts)
    return(matrix)
    
  }else{
    rows <- NROW(matrix)
    diversity <- c()
    dispersion <- c()
    restart <- c()

    for(i in 1:rows){
      current_div <- get_diversity(matrix[i,], distances)
      diversity <- c(diversity, current_div)
      current_dis <- get_dispersion(matrix[i,], distances)
      dispersion <- c(dispersion, current_dis)
      restart <- c(restart, restarts)
    }
    matrix <- cbind(diversity, dispersion, restarts)
    return(matrix)
  }
}


get_diversity <- function(partition, distances){
  
  sum = 0
  n = length(partition)
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(partition[i] == partition[j]){
        sum = sum + distances[i,j]
      }
    }
  }
  
  return(sum);
}


get_dispersion <- function(partition, distances){
  
  min = Inf
  n = length(partition)
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(partition[i] == partition[j]){
        distance = distances[i,j];
        if(distance < min){
          min = distance;
        }
      }
    }
  }
  
  return(min);
}