#converts solution of ILP into a partition
#########################
  assignement_to_partition <- function(assignement){
    assignement <- assignement[,-c(1,4)]
    partition <- assignement[order(assignement$i),]
    partition <- partition[,-c(1)]
    return(partition)
  }
#create random matrix with 3 attributes
############################################
row1 <- matrix(runif(6,0,100),6,1)
row2 <- matrix(runif(6,0,50),6,1)
row3 <- matrix(runif(6,0,30),6,1)
matrix <- cbind(row1,row2,row3)
############################################
#cres2 <- bicriterion_iterated_local_search_call(matrix, 2,2)
#ilpsol <- cluster_assignement_ilp(matrix,2)

#create multiple exchange solutions
#######################################
# 
#exchange_solution <- anticlustering(
#rmatrix,
#K = 10,
#objective = "distance",   method = 
#)
#
# exchange_solutions <- c()
# exchange_solutions <- rbind(exchange_solutions,exchange_solution)
# exchange_plot <- plot_partition_matrix(exchange_solutions, rmatrix, 0)
# View(exchange_plot)
#  
# start_time <- Sys.time()
#cres2 <- bicriterion_iterated_local_search_call(matrix, 10,2)
#plot2 <- plot_partition_matrix(cres2, matrix, 2)
# both_heuristics <- rbind(exchange_plot, plot2)
# both_heuristics <- as.data.frame(both_heuristics)
#  
# ggplot(both_heuristics, aes(x = diversity, y = dispersion, color = factor(restarts, labels = c("exchange method", 
#  "2-restart-BILS")))) + geom_point() + labs(color = "heuristic:")
#
  
#plot different restarts of the BILS heuristic
###################################################################  
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
# combi_plot <- rbind(plot1, plot10, plot100, plot1000, plot10000)
#
#  
# combi_plot2 <- as.data.frame(combi_plot2)
# ggplot(combi_plot2, aes(x = diversity, y = dispersion, shape = factor(restarts), color = factor(restarts))) + 
# geom_point() + labs(shape = "restarts", color = "restarts")
#
# ggplot(combi_plot, aes(x = diversity, y = dispersion, shape = factor(restarts), color = factor(restarts))) + 
# geom_point() + labs(shape = "restarts", color = "restarts")
#
# ggplot(combi_plot, aes(x = diversity, y = dispersion)) + geom_point(aes(color = factor(restarts))) + theme(aspect.ratio=1) + 
# labs(color = "Number of restarts")
# 
#
  
# check for partitions with the same diversity
########################################################
# 
# N <- 10
# K <- 2
# features <- matrix(sample(N, replace = TRUE)) # Ganzzahlige Zufallsdaten
# features <- matrix(rnorm(N))
# partitions <- generate_partitions(N, K)
# length(partitions) # number of possible partitions
# 
# Create an objective function that takes the partition
# as first argument (then, we can use sapply to compute
# the objective for each partition)
# var_obj <- function(clusters, features) {
#   distance_objective(features, clusters)
# } 
# 
#   all_objectives <- sapply(
#     partitions,
#     FUN = var_obj,
#     features = features
#   )
# 
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