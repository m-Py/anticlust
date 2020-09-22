library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)
library(knitr)


cluster_assignement_ilp <- function(data, g){

  distances <- convert_to_distances(data)
  n <- NROW(distances)
  m <- max(distances)


  model <- MIPModel() %>%

    add_variable(x[i, k], i = 1:n, k = 1:g, type = "binary") %>%

    add_variable(z[i, j, k], i = 1:n, j = 1:n, k = 1:g, type = "binary") %>%

    add_variable(min_distance, type = "continuous") %>%

    set_objective(min_distance, sense = "max") %>%

    add_constraint(sum_expr(x[i, k], k = 1:g) == 1, i = 1:n) %>%

    add_constraint(min_distance <= distances[i,j] * z[i,j,k] + (1 - z[i, j, k]) * m , i < j, i = 1:n, j = 1:n, k = 1:g) %>%

    add_constraint(x[i, k] + x[j, k] <= 1 + z[i, j ,k], i < j, i = 1:n, j = 1:n, k = 1:g) %>%
    
    add_constraint(sum_expr(x[i, k], i = 1:n) == n/g, k = 1:g) %>%

    set_bounds(z[i, j, k], i >= j, ub = 0, i = 1:n, j = 1:n, k = 1:g)
  model
  View(model)
  result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
  solution <- get_solution(result, x[i,k]) %>%
    filter(value > 0)

  return (solution)
}