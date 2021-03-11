#' Solve anticluster editing using MBPI
#'
#' @param data An N x N dissimilarity matrix
#' @param G The number of clusters 
#' @param R The desired number of restarts for the algorithm
#' @param W list of possible weights for dispersion in the bicriterion (0 <= W <= 1)
#'
#' @return pareto set
#' @export


multistart_bicriterion_pairwise_interchange <- function(data, G, R, W){ 
  distances <- convert_to_distances(data) 
  pareto_set = list()
  for(i in 1:R){
    w1 <- sample(W, 1)
    w2 <- 1 - w1
    partition <- sample(rep_len(1:G, length.out = nrow(distances)))
    pareto_set <- update_pareto(partition, pareto_set, distances)
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
  }
  return(pareto_set)
}

sum_distances <- function(partition, distances){
  sum <- 0
  for(i in 1:(length(partition)-1)){
    for(j in (i+1):length(partition)){
      if(partition[i] == partition[j]){
        sum <- sum + distances[i,j]
      }
    }
  }
  return(sum)
}

minimal_distance <- function(partition, distances){
  min <- Inf
  for(i in 1:(length(partition)-1)){
    for(j in (i+1):length(partition)){
      if(partition[i] == partition[j]){
        distance <- distances[i,j]
        if(distance < min){
          min <- distance
        }
      }
    }
  }
  return(min)
}


update_pareto <- function(partition, pareto_set, data){
  diversity <- sum_distances(partition, data)
  dispersion <- minimal_distance(partition, data)
  partition <- compress(partition)
  pareto_element <- list(partition, diversity, dispersion)
  
  # in the first run the pareto_set is empty and the element is added directly
  if(length(pareto_set) == 0){
    pareto_set[[length(pareto_set)+1]] <- pareto_element
    return(pareto_set)
  }
  
  dominated <- paretodominated(diversity, dispersion, pareto_set)
  included <- pareto_included(pareto_element, pareto_set)
  if(!dominated & !included){
    pareto_set <- delete_outdated(diversity, dispersion, pareto_set)
    pareto_set[[length(pareto_set)+1]] <- pareto_element
    
  }  
  return(pareto_set)
}

#a partition is dominated by the incumbent set
#if one elementof it is better in both diversity and dispersion
paretodominated <- function(diversity, dispersion, pareto_set){
  for(i in 1:length(pareto_set)){
    if((pareto_set[[i]][[2]] >= diversity & pareto_set[[i]][[3]] > dispersion) | (pareto_set[[i]][[2]] > diversity & pareto_set[[i]][[3]] >= dispersion)){
      return(TRUE)
    }
  } 
  return (FALSE)
} 

#this exact partition is already in the set
pareto_included <- function(pareto_element, pareto_set){
  for(entry in pareto_set){
    if(all(entry[[1]] == pareto_element[[1]])){
      return(TRUE)
    }
  }
  return(FALSE)
}

#elements in the old paretoset dominated by a new one are deleted from the set
delete_outdated <- function(diversity, dispersion, pareto_set){
  i <- 1
  while(i <= length(pareto_set)){
    if((diversity >= pareto_set[[i]][[2]] & dispersion > pareto_set[[i]][[3]]) | (diversity > pareto_set[[i]][[2]] & dispersion >= pareto_set[[i]][[3]])){
      pareto_set[[i]] <- NULL
      i <- i - 1
    }
    i <- i + 1
  }
  return (pareto_set)
}

#partitions with the same clusters but different numbers
#for them get compressed into the same ascending format
compress <- function(vector){
  cluster <- 1
  for(i in 1:length(vector)) {
    if(vector[i] < cluster){
      next
    }else if(vector[i] == cluster){
      cluster <- cluster + 1
      next
    }
    change <- vector[i]
    if(cluster < change){
      for(j in 1:length(vector)){
        if(change == vector[j]){
          vector[j] <- cluster
        }else if(cluster == vector[j]){
          vector[j] <- change
        }
      }
    }
    cluster <- cluster + 1
  }
  return(vector)
}