library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)
library(knitr)


dispersion_ilp <- function(data, g){
  
  distances <- convert_to_distances(data)
  n <- NROW(distances)
  m <- max(distances)
  
  
  model <- MIPModel() %>%
    
    add_variable(x[i, j], i = 1:n, j = 1:n, type = "binary") %>%
    
    add_variable(min_distance, type = "continuous") %>%
    
    set_objective(min_distance, sense = "max") %>% 
    
    add_constraint(min_distance <= distances[i,j] * x[i,j] + (1 - x[i,j]) * m , i = 1:n, j = 1:n) %>%
    
    add_constraint(x[i, k] + x[j, k] - x[i, j] <= 1, i = 1:n, j = 1:n, k = 1:n, i < j, i < k, j < k) %>%
    
    add_constraint(x[i, j] + x[j, k] - x[i, k] <= 1, i = 1:n, j = 1:n, k = 1:n, i < j, i < k, j < k) %>%
    
    add_constraint(x[i, j] + x[i, k] - x[j, k] <= 1, i = 1:n, j = 1:n, k = 1:n, i < j, i < k, j < k) %>%
    
    add_constraint(sum_expr(x[i, j], i < j, j = 1:n) + sum_expr(x[k, i], k < i, k = 1:n) == n/g - 1, i = 1:n) %>% 
    
    set_bounds(x[i, j], i >= j, ub = 0, i = 1:n, j = 1:n)
  model
  
  result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
  solution <- get_solution(result, x[i, j]) %>%
    filter(value > 0)
  
  return (solution)
}



