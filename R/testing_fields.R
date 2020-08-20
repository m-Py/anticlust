  convert_partition <- function(partition){
    partition_list <- list()
    for(i in 1:max(partition)){
      cluster <- c()
      pos <- 1
      for(j in partition){
        if(i == j){
          cluster <- c(cluster, pos)
        }
        pos <- pos + 1
      }
      partition_list[[length(partition_list)+1]] <- cluster
    }
    return ( partition_list)
  }
  
  vectorize <- function(list, index){
    vector <- c()
    for(element in list){
      vector <- c(vector, element[[index]])
    }
    return(vector)
  }
  
  plot_set <- function(pareto_set){
    diversity <- vectorize(pareto_set, 2)
    dispersion <- vectorize(pareto_set, 3)
    plot(diversity, dispersion)
  }
  
  
  add_element <- function(partition, pareto_set, data){
    diversity <- sum_distances(partition, convert_to_distances(data))
    dispersion <- minimal_distance(partition, convert_to_distances(data))
    partition <- compress(partition)
    pareto_element <- list(partition, diversity, dispersion)
    pareto_set[[length(pareto_set)+1]] <- pareto_element
    return(pareto_set)
  }
  
  edges_to_partition <- function(edges){
    edges <- edges[,-c(1,4)]
    partition <- vector(length =max(edges))
    counter <- 1
    while(NROW(edges) > 0){
      min <- min(edges)
      pos <- 1
      while(pos <= NROW(edges)){
        if(edges[[1]][[pos]] == min | edges[[2]][[pos]] == min | partition[edges[pos,1]] == counter | partition[edges[pos,2]] == counter){
          partition[edges[pos,1]] <- counter
          partition[edges[pos,2]] <- counter
          edges <- edges[-pos,]
          pos <- 1
          next
        }
        pos <- pos + 1
      }
      counter <- counter + 1
    }
    return (partition)
  }
  
  assignement_to_partition <- function(assignement){
    assignement <- assignement[,-c(1,4)]
    partition <- assignement[order(assignement$i),]
    partition <- partition[,-c(1)]
    return(partition)
  }
  
  
    

matrix <- rbind(c(0,1,2,1,3,4,3,3,1),
                c(1,0,2,3,4,1,2,2,3),
                c(2,2,0,1,1,4,3,2,1),
                c(1,3,1,0,3,2,1,3,4),
                c(3,4,1,3,0,1,1,2,1),
                c(4,1,4,2,1,0,3,2,3),
                c(3,2,3,1,1,3,0,1,2),
                c(3,2,2,3,2,2,1,0,4),
                c(1,3,1,4,1,3,2,4,0))

matrix2 <- rbind(c(0,1,2,1,3,4),
                 c(1,0,2,3,4,1),
                 c(2,2,0,1,1,4),
                 c(1,3,1,0,3,2),
                 c(3,4,1,3,0,1),
                 c(4,1,4,2,1,0))

matrix3 <- rbind(c(0,1,2,3),
                 c(1,0,2,3),
                 c(2,2,0,1),
                 c(3,3,1,0))

  
#x <- diversity_ilp(matrix3, 2)

  # dumdum <- bicriterion_iterated_local_search_call(mymatrix, 2, 100)
  # plot_partition_matrix(dumdum,mymatrix,100)
  # View(dumdum)
  # 
  #     # # combi_result <- add_element(div2_partition, r_result, matrix2)
  # # combi_result <- add_element(div2_partition, combi_result, matrix)
  # # View(combi_result)
  # # # # # combi_result <- add_element(disp_partition,combi_result, matrix)
  # # # View(combi_result)
  # # # plot_set(combi_result)
  # # # # 
  # # # 
 #curres <- cluster_assignement_ilp(matrix3, 2)
  # # #View(curres)
  # # curres <- assignement_to_partition(curres)
  # # together <- add_element(curres, r_result, matrix)
  # # View(together)
  # #View(partition)
  # 
  # 
  # #ilp <- add_element(parti, lowkey  , matrix)
  # #plot_set(r_result)
  # #plot_set(ilp)
  # #View(combi_result)
  # # plot_partition_matrix(c_result, matrix2)
  # 
  # 
  # 
  #
# library(ggplot2)
row1 <- matrix(runif(50,0,100),50,1)
row2 <- matrix(runif(50,0,50),50,1)
row3 <- matrix(runif(50,0,30),50,1)
matrix <- cbind(row1,row2,row3)
# View(matrix)
  # rmatrix <- convert_to_distances(rmatrix)
  # #View(rmatrix)
  # exchangemethod <- function(){
  #   start_time <- Sys.time()
  #   for(i in (1:10)){
  #     exchange_solution <- anticlustering(
  #     rmatrix,
  #     K = 10,
  #     objective = "distance",   method = "exchange"
  #     )
  #   }
  #   end_time <- Sys.time()
  #   print(end_time - start_time)
  #   return (exchange_solution)
  #   }
  # x = exchangemethod()
  # View(x)
  # 
  # 
#exchange_solution <- anticlustering(
#rmatrix,
#K = 10,
#objective = "distance",   method = 
#)

#x <- anticlustering(matrix, 5, objective = "distance", method ="ilp")
#y <- bicriterion_iterated_local_search_call(matrix, 3, 10)
# View(x)
# View(y)
# plot_partition_matrix(y, matrix, 1)
#   # 
  # #exchange_solutions <- x
  # exchange_solutions <- rbind(exchange_solutions,exchange_solution)
  # exchange_plot <- plot_partition_matrix(exchange_solutions, rmatrix, 0)
  # View(exchange_plot)
  #  
  # start_time <- Sys.time()
#cres300000 <- bicriterion_iterated_local_search_call(matrix, 10,300000)
cres10000 <- bicriterion_iterated_local_search_call(matrix, 10,10000)
  plot300000 <- plot_partition_matrix(cres300000, matrix, 300000)
  plot10000 <- plot_partition_matrix(cres10000, matrix, 300000)
  View(plot300000)
  View(plot10000)
  # print(end_time - start_time)
  #     plot1 <- plot_partition_matrix(cres1, rmatrix, 10000)
  # View(plot1)
  # #View(plot2)
  # # # # cres16 <- bicriterion_iterated_local_search_call(rmatrix, 10, 16)
  # # # # cres10000 <- bicriterion_iterated_local_search_call(rmatrix, 10, 10000)
  # # plot1 <- plot_partition_matrix(cres1, rmatrix, 1)
  # both_heuristics <- rbind(exchange_plot, plot1)
  # both_heuristics <- as.data.frame(both_heuristics)
  # # # 
  # ggplot(both_heuristics, aes(x = diversity, y = dispersion, color = factor(restarts, labels = c("exchange method", 
  #  "2-restart-BILS")))) + geom_point() + labs(color = "heuristic:")
  # 
    # cres1 <- bicriterion_iterated_local_search_call(rmatrix, 10, 1)
  # cres10 <- bicriterion_iterated_local_search_call(rmatrix, 10, 10)
  # cres100 <- bicriterion_iterated_local_search_call(rmatrix, 10, 100)
  # cres1000 <- bicriterion_iterated_local_search_call(rmatrix, 10, 1000)
  # cres10000 <- bicriterion_iterated_local_search_call(rmatrix, 10, 10000)
  # plot1 <- plot_partition_matrix(cres1, rmatrix, 1)
  # plot10 <- plot_partition_matrix(cres10, rmatrix, 10)
  # plot100 <- plot_partition_matrix(cres100, rmatrix, 100)
  # plot1000 <- plot_partition_matrix(cres1000, rmatrix, 1000)
  # plot10000 <- plot_partition_matrix(cres10000, rmatrix, 10000)
  # combi_plot2 <- rbind(plot1, plot10, plot100, plot1000, plot10000)
  # # 
  # #  
  # combi_plot2 <- as.data.frame(combi_plot2)
    # ggplot(combi_plot2, aes(x = diversity, y = dispersion, shape = factor(restarts), color = factor(restarts))) + 
    #   geom_point() + labs(shape = "restarts", color = "restarts")
  # ggplot(combi_plot, aes(x = diversity, y = dispersion, shape = factor(restarts), color = factor(restarts))) + 
  # geom_point() + labs(shape = "restarts", color = "restarts")
  #   # ggplot(combi_plot, aes(x = diversity, y = dispersion)) + geom_point(aes(color = factor(restarts))) + theme(aspect.ratio=1) + 
  # #   labs(color = "Number of restarts")
  # #  
  # combi_plot <- rbind(plot1, plot10, plot100, plot1000, plot10000)
  #  
  # 
#    ilprow1 <- matrix(runif(20,0,100),20,1)
#    ilprow2 <- matrix(runif(20,0,50),20,1)
#    ilprow3 <- matrix(runif(20,0,10),20,1)
#    ilpmatrix <- cbind(ilprow1, ilprow2, ilprow3)
#  View(ilpmatrix)
#disp_ilp1 <- dispersion_ilp(matrix, 3)
# disp_partition <- edges_to_partition(disp_ilp1)
#disp_ilp2 <- cluster_assignement_ilp(matrix, 3)
#   # disp_partition2 <- assignement_to_partition(disp_ilp2)
# View(disp_ilp1)
  # compare <- c()
  # compare <- add_element(disp_partition, compare, ilpmatrix)
  # compare <- add_element(disp_partition2, compare, ilpmatrix )
  #^ View(compare)
  # res <- bicriterion_iterated_local_search_call(ilpmatrix,2, 100)
  # plot_partition_matrix(res, ilpmatrix,5)
  # plot_set(compare)
  # View(compare)
  # 
  # # ######################################################
  # 
  # N <- 10
  # K <- 2
  # features <- matrix(sample(N, replace = TRUE)) # Ganzzahlige Zufallsdaten
  # features <- matrix(rnorm(N))
  # partitions <- generate_partitions(N, K)
  # length(partitions) # number of possible partitions
  # 
  # ## Create an objective function that takes the partition
  # ## as first argument (then, we can use sapply to compute
  # ## the objective for each partition)
  # var_obj <- function(clusters, features) {
  #   distance_objective(features, clusters)
  # }
  # 

# 
#   all_objectives <- sapply(
#     partitions,
#     FUN = var_obj,
#     features = features
#   )
# 
#   x <-round(all_objectives,3)
#   sum(x == max(x))
#   all_objectives[max(all_objectives)]
#   View(all_objectives)
# 
#   unique(all_objectives)
# 
#   # # Der Anteil der Partitionen, die dasselbe Objective haben
#   length(all_objectives)
#   mean(all_objectives == max(all_objectives))
  ############################################################
  
  # pareto_set2 <- function(df) {
  #   df <- data.frame(df)
  #   ncrits <- ncol(df)
  #   df$ID <- 1:nrow(df)
  #   df <- df[order(rowSums(df[, 1:ncrits]), decreasing = TRUE), ]
  #   
  #   is_dominated <- rep(FALSE, nrow(df)) 
  #   for (i in nrow(df):1) {
  #     for (j in 1:nrow(df)) {
  #       j_dominates_i <- all(df[i, 1:ncrits] <= df[j, 1:ncrits]) && 
  #         any(df[i, 1:ncrits] < df[j, 1:ncrits])
  #       if (j_dominates_i) {
  #         is_dominated[df$ID[i]] <- TRUE
  #         break
  #       }
  #     }
  #   }
  #   !is_dominated
  # }
  # 
  # library(anticlust)
  # 
  # N <- 9 # number of elements
  # K <- 3 # number of groups
  # M <- 4 # number of variables
  # all_partitions <- generate_partitions(N, K) # all partitions
  # features <- matrix(rnorm(N * M), ncol = M) # data
  # View(all_partitions)
  # 
  # 
  # diversity <- function(partitions, data) {
  #   anticlust::diversity_objective(data, partitions)
  # }
  # 
  # dispersion <- function(partitions, data) {
  #   anticlust::dispersion_objective(data, partitions)
  # }
  # 
  # # For each partition, compute the diversity objective
  # all_objectives_diversity <- sapply(
  #   all_partitions,
  #   FUN = diversity,
  #   data = features
  # )
  # 
  # # For each partition, compute the dispersion objective
  # all_objectives_dispersion <- sapply(
  #   all_partitions,
  #   FUN = dispersion,
  #   data = features
  # )
  # 
  # # Store objectives as data frame
  # df <- data.frame(
  #   diversity = all_objectives_diversity,
  #   dispersion  = all_objectives_dispersion
  # )
  # 
  # # Now: Compute Pareto efficient set
  # 
  # # For a data frame (rows = partitions, columns = objectives), determine which
  # # elements are part of the Pareto efficient set. Output is a logical vector,
  # # illustrating for each row if it is part of the Pareto efficient set.
  # # (function is very slow / inefficient)
  # pareto_set <- function(df) {
  #   is_dominated <- rep(FALSE, nrow(df))
  #   for (i in 1:nrow(df)) {
  #     for (j in 1:nrow(df)) {
  #       j_dominates_i <- all(df[i, ] <= df[j, ]) && any(df[i, ] < df[j, ])
  #       if (j_dominates_i) {
  #         is_dominated[i] <- TRUE
  #         break
  #       }
  #     }
  #   }
  #   !is_dominated
  # }
  # 
  # efficient_set <- pareto_set(df)
  # 
  # # Plot Pareto front
  # plot(df, col = "darkgrey", las = 1, pch = 4, cex = 0.8)
  # points(df[efficient_set, ], col = "red", cex = 1.2, pch = 19)
  # 
