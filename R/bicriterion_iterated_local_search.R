#' Solve anticluster editing using BILS
#'
#' @param data An N x N dissimilarity matrix
#' @param G The number of clusters 
#' @param R The desired number of restarts for the algorithm
#' @param W list of possible weights for dispersion in the bicriterion (0 <= W <= 1)
#' @param Xi a vector for the neighborhood size range [xi1,xi2] in percent 
#' (for example c(0.05, 0.1) can be interpreted as 5%-10%)
#'
#' @return pareto set
#' @export

bicriterion_iterated_local_search <- function(data, G, R, W, Xi){
  distances <- convert_to_distances(data) 
  pareto_set <- multistart_bicriterion_pairwise_interchange(data, G, (ceiling(R/2)), W)
  for(q in 1:(ceiling(R/2))){
    w1 <- sample(W, 1)
    w2 <- 1 - w1
    neighborhood_size <- runif(1, Xi[1], Xi[2])
    pareto_sample <- sample(pareto_set,1)
    partition <- pareto_sample[[1]][[1]]
    for(i in 1:(nrow(distances)-1)){
      for(j in (i+1):(nrow(distances))){
        g <- partition[i]
        h <- partition[j]
        if(g != h){
          random <- runif(1, 0, 1)
          if(random < neighborhood_size){
            partition <- cluster_swap(partition,i,j)
          }
        }
      }
    }
    pareto_set <- local_search(partition, distances, pareto_set, w1, w2)
  }
  return(pareto_set)
} 


local_search <- function(partition, distances, pareto_set, w1, w2){
  diversity <- sum_distances(partition, distances)
  dispersion <- minimal_distance(partition, distances)
  max_bicriterion <- w1*diversity + w2*dispersion
  Flag <- FALSE
  while(!Flag){
    Flag <- TRUE
    for(i in 1:(nrow(distances)-1)){
      for(j in (i+1):nrow(distances)){
        g <- partition[i]
        h <- partition[j]
        if(g != h){
          partition <- cluster_swap(partition,i,j)
          pareto_set <- update_pareto(partition, pareto_set, distances)
          current_diversion <- sum_distances(partition, distances)
          current_dispersion <- minimal_distance(partition, distances)
          current_bicriterion <- w1*current_diversion + w2*current_dispersion
          if(current_bicriterion > max_bicriterion){
            max_bicriterion <- current_bicriterion
            Flag <- FALSE
          }else{
            partition <- cluster_swap(partition,i,j)
          }
        }
      }
    }
  }
  return(pareto_set)
}