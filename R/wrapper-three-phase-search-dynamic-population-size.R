#' Three phase search with dynamic population size heuristic
#'
#' This function implements the three phase search algorithm TPSPD for
#' anticlustering by Yang et al. (2022; <doi.org/10.1016/j.ejor.2022.02.003>).
#' The description of their algorithm is
#' given in Section 2 of their paper (in particular, see the
#' Pseudocode in Algorithm 1).
#' 
#' 
#' @param x The data input, as in \code{\link{anticlustering}}.
#' @param K Number of anticlusters to be formed.
#' @param N Number of elements.
#' @param objective The anticlustering objective, can be "diversity" or "dispersion".
#' @param number_iterations A number that defines how many times the steps in the search algorithm are repeated.
#' @param clusters A vector of length K that specifies the number of elements each cluster can contain. 
#' If this vector is not NULL, the lower and upper bounds will be disregarded.
#' @param beta_max The algorithm begins with a pool of random initial solutions of size beta_max. 
#'  Over time, the size of the solution pool decreases linearly until it reaches beta_min.
#' @param beta_min The minimum solution pool size the algorithm should reach before making a determination.
#' @param lower_bound Minimum number of elements in each anticluster. By default, anticlusters are of equal size,
#'  calculated as the total number of items divided by the number of clusters.
#' @param upper_bound Maximum number of elements in each anticluster. By default, anticlusters are of equal size,
#'  calculated as the total number of items divided by the number of clusters.
#' @param theta_max Parameter for the strength of undirected perturbation, 
#' which decreases linearly over time from theta_max to theta_min.
#' @param theta_min Parameter for the strength of undirected perturbation, 
#' which decreases linearly over time from theta_max to theta_min.
#' @param eta_max Parameter that specifies how many times the steps in the direct perturbation are executed.
#' @param alpha Parameter for weighing the discrimination of a slightly worse local optimal child solution.
#'     
#' @details Details of the implementation of the algorithm can be found 
#'  in the pseudocode of the paper Yang et al. (2022). However, we performed one change
#'  as compared to the original description of the algorithm: Instead of 
#'  setting a time limit, we define the number of iterations the algorithm 
#'  performs before terminating (via argument \code{number_iterations}).
#' 
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and \code{K}) to each input element
#' 
#' @author Hannah Hengelbrock \email{Hannah.Hengelbrock@@hhu.de}, 
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @export
#' 
#' @examples 
#' 
#' # Generate some random data
#' N <- 120
#' M <- 5
#' K <- 4
#' dat <- matrix(rnorm(N * M), ncol = M)
#' distances <- dist(dat)
#'
#' # Perform three hase serach algorithm
#' result1 <- three_phase_search_anticlustering(dat, K, N)
#'
#' # Compute objectives function
#' diversity_objective(distances, result1)
#' 
#' # Standard algorithm:
#' result2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
#' diversity_objective(distances, result2)
#' 
#' 
#' @references
#' 
#' Xiao Yang et al. “A three-phase search approach with dynamic population size for solving 
#' the maximally diverse grouping problem”. In: European Journal of Operational Research
#' 302.3 (2022) <doi:10.1016/j.ejor.2022.02.003>
#' 
three_phase_search_anticlustering <- function(x, K, N, objective = "diversity",
    number_iterations=50, clusters = NULL, upper_bound  = NULL, lower_bound  = NULL, 
    beta_max = 15,  theta_max = NULL, theta_min = NULL, beta_min = NULL, eta_max=3, alpha=0.05) {

    input_validation_threephase_search(x, K, N, objective, clusters, number_iterations, upper_bound, 
    theta_max, theta_min, lower_bound, beta_max, beta_min, eta_max, alpha)
    
    distances <- convert_to_distances(x) 

    if (is.null(lower_bound)) {
       lower_bound <- floor(N/K)
    } 
    if (is.null(upper_bound)) {
       upper_bound <- ceiling(N/K)
    }   
    
    if (N <= 400  & is.null(theta_max) & is.null(theta_min) & is.null(beta_min)) {
    	theta_max  <- 1.2
    	theta_min  <- 0.1
    	beta_min  <- 2
    } else if ( is.null(theta_max) & is.null(theta_min) & is.null(beta_min)) {
    	theta_max  <- 2.0
    	theta_min  <- 1.0
    	beta_min  <- 1
    }

     # create result vector for results to use in C
    result_vector <- numeric(N)
    
    # If clusters are not predefined, set the value to -1 so the C implementation knows 
    # to evenly distribute the number of elements based on K and the boundaries.
    if (is.null(clusters)) {
      clusters <- rep(-1, K)
    }
     
     if (objective == "diversity") {
     results <- .C("three_phase_search_dynamic_population_size",
                  distances = as.double(distances),
                  N_in = as.integer(N),
                  K_in = as.integer(K),
                  number_of_iterations = as.integer(number_iterations),
                  clusters = as.integer(clusters),
                  upper_bound = as.integer(upper_bound),
                  lower_bound = as.integer(lower_bound),
                  Beta_max = as.integer(beta_max),
                  elapsed_time = as.integer(0),
                  Theta_max = as.double(theta_max),
                  Theta_min = as.double(theta_min),
                  Beta_min = as.integer(beta_min),
                  Eta_max = as.integer(eta_max),
                  Alpha = as.double(alpha),
                  result = as.integer(result_vector),
                  score = as.double(0.0),
                  mem_error = as.integer(0),
                  PACKAGE = "anticlust"
                  )
     } else if (objective == "dispersion") {
      results <- .C("three_phase_search_dispersion",
                  distances = as.double(distances),
                  N_in = as.integer(N),
                  K_in = as.integer(K),
                  number_of_iterations = as.integer(number_iterations),
                  clusters = as.integer(clusters),
                  upper_bound = as.integer(upper_bound),
                  lower_bound = as.integer(lower_bound),
                  Beta_max = as.integer(beta_max),
                  elapsed_time = as.integer(0),
                  Theta_max = as.double(theta_max),
                  Theta_min = as.double(theta_min),
                  Beta_min = as.integer(beta_min),
                  Eta_max = as.integer(eta_max),
                  Alpha = as.double(alpha),
                  result = as.integer(result_vector),
                  score = as.double(0.0),
                  mem_error = as.integer(0),
                  PACKAGE = "anticlust"
                  )
     } else {
      stop("Objective funtion for TPSDP is not defined.")
     }

     results[["mem_error"]]
     if (results[["mem_error"]] == 1) {
       stop("Could not allocate enough memory.")
     }

    return(results$result + 1)
}

input_validation_threephase_search <- function(x, K, N, objective, clusters, number_iterations, upper_bound, lower_bound, 
theta_max, theta_min, beta_max, beta_min, eta_max, alpha) {

    # cluster vector
    if (!is.null(clusters)) {
    validate_input(clusters, "clusters", len = K)
    }
      # Objective
    validate_input(
    objective, "objective", objmode = "character", len = 1,
    input_set = c("diversity", "dispersion"), not_na = TRUE,  not_function = TRUE
    )
    validate_input(K, "K",  must_be_integer = TRUE, not_na = TRUE)
    validate_input(lower_bound, "lower_bound", greater_than = 0, must_be_integer = TRUE)
    validate_input(upper_bound, "upper_bound", greater_than = 0, must_be_integer = TRUE)
    validate_input(beta_max, "beta_max", greater_than = 0, must_be_integer = TRUE)
    validate_input(beta_min, "beta_min", greater_than = 0, must_be_integer = TRUE)
    validate_input(eta_max, "eta_max", greater_than = 0, must_be_integer = TRUE)
    validate_input(theta_max, "theta_max", greater_than = 0.0)
    validate_input(theta_min, "theta_min", greater_than = 0.0)
    validate_input(eta_max, "eta_max", greater_than = 0, must_be_integer = TRUE)
    validate_input(alpha, "alpha", greater_than = 0.0)
}
