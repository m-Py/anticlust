#' Feasible and Infeasible Region search heuristic
#'
#' This function implements the feasible and infeasible region search algorithm FIFR for
#' anticlustering by Wu et al. (2025; <doi.org/10.1016/j.cor.2025.107030>).
#' The description of their algorithm is given in Section 4 of their paper (in particular, see the
#' Pseudocode in Algorithm 1).
#' 
#' 
#' @param matrix The data input. Currently just a vector.
#' @param K Number of anticlusters to be formed.
#' @param N Number of elememts.
#' @param objective The anticlustering objective = "diversity".
#' @param number_iterations A number that defines how many times the steps in the search algorithm are repeated.
#' @param clusters A vector of length M that specifies the number of elements each cluster can contain. 
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
#' @param phi Parameter that determines the population size when initiating a new round of exploration upon
#' triggering the jump-back mechanism
#' @param tau Parameter when "noImp" exceeds it, the algorithm transits from exploitation to exploration strategy
#' @param kmax A parameter that determines the maximum degree of constraint violation in the IFR search
#' @param alpha Parameter for weitghing the discrimitation of a slighlty worse local optiomal child solution
#'     
#' @details Details of the implementation of the algorithm can be found 
#'  in the pseudocode of the paper Wu et al. (2025)
#' 
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and \code{K}) to each input element.
#' 
#' @author David Buczynski \email{david.buczynski@@hhu.de}, 
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @export
#' 
#' @examples 
#' 
#' # Generate some random data
#' N <- 120
#' M <- 5
#' K <- 3
#' dat <- matrix(rnorm(N * M), ncol = M)
#' distances <- dist(dat)
#'
#' # Perform three hase serach algorithm
#' results1 <- feasible_and_infeasible_region_search_anticlustering(distances, K, N)
#' results2 <- anticlustering(distances, K = K, method = "3phase")
#' results3 <- anticlustering(distances, K = K, method = "local-maximum", repetitions = 50)
#'
#' # Compute objectives funtion
#' diversity_objective(distances, results1)
#' diversity_objective(distances, results2)
#' diversity_objective(distances, results3)
#' 
#' 
#' @references
#' 
#' Wu, X., Feng, J., Yang, J., & Zhang, Y. (2025). Feasible and infeasible region 
#' search for the maximally diverse grouping problem. Computers & Operations Research, 179,
#' 107030. https://doi.org/10.1016/j.cor.2025.107030
#'
feasible_and_infeasible_region_search_anticlustering <- function(
    x, K, N, number_iterations=50, clusters=NULL, upper_bound=NULL, 
    lower_bound=NULL, beta_max=NULL, theta_max=NULL, theta_min=NULL, beta_min=NULL, phi=0.8,
    tau=NULL, kmax=4, alpha=NULL
) {

    distances <- convert_to_distances(x)

    if (is.null(lower_bound)) {
        lower_bound <- floor(N/K)
    }
    if (is.null(upper_bound)) {
        upper_bound <- ceiling(N/K)
    }

if (N < 480) {
    if (is.null(theta_max)) theta_max <- 1.2
    if (is.null(theta_min)) theta_min <- 0.1
    if (is.null(alpha)) alpha <- 0.1
    if (is.null(tau)) tau <- 300
    if (is.null(beta_max)) beta_max <- 10
    if (is.null(beta_min)) beta_min <- 1
} else if (N < 960) {
    if (is.null(theta_max)) theta_max <- 2.0
    if (is.null(theta_min)) theta_min <- 0.1
    if (is.null(alpha)) alpha <- 0.05
    if (is.null(tau)) tau <- 20
    if (is.null(beta_max)) beta_max <- 10
    if (is.null(beta_min)) beta_min <- 1
} else {
    if (is.null(theta_max)) theta_max <- 2.0
    if (is.null(theta_min)) theta_min <- 0.1   
    if (is.null(alpha)) alpha <- 0.05          
    if (is.null(tau)) tau <- 20
    if (is.null(beta_max)) beta_max <- 10
    if (is.null(beta_min)) beta_min <- 0
}


    input_validation_fifr_search(x, K, N, clusters, number_iterations, upper_bound, 
    lower_bound, theta_max, theta_min, beta_max, beta_min, phi, tau, kmax, alpha)

    # create result vector for results to use in C
    result_vector <- numeric(N)

    # If clusters are not predefined, set the value to -1 so the C implementation knows 
    # to evenly distribute the number of elements based on K and the boundaries.
    if (is.null(clusters)) {
      clusters <- rep(-1, K)
    }

   results <- .C("feasible_and_infeasible_region_search",
                  distances = as.double(distances),
                  N_in = as.integer(N),
                  M_in = as.integer(K),
                  number_of_iterations = as.integer(number_iterations),
                  clusters = as.integer(clusters),
                  upper_bound = as.integer(upper_bound),
                  lower_bound = as.integer(lower_bound),
                  Beta_max = as.integer(beta_max),
                  elapsed_time = as.integer(0),
                  Theta_max = as.double(theta_max),
                  Theta_min = as.double(theta_min),
                  Beta_min = as.integer(beta_min),
                  Phi = as.double(phi),
                  Tau = as.integer(tau),
                  Kmax = as.integer(kmax),
                  Alpha = as.double(alpha),
                  result = as.integer(result_vector),
                  score = as.double(0.0),
                  mem_error = as.integer(0),
                  PACKAGE = "anticlust"
    )
    
    results[["mem_error"]]
    if (results[["mem_error"]] == 1) {
       stop("Could not allocate enough memory.")
    }
    
    results$result + 1 # in C, we use 0, 1, 2... for cluster labels
}

input_validation_fifr_search <- function(
    x, K, N, clusters, number_iterations, upper_bound, lower_bound, 
    theta_max, theta_min, beta_max, beta_min, phi, tau, kmax, alpha
) {

    # cluster vector
    if (!is.null(clusters)) {
        validate_input(clusters, "clusters", len = K)
    }

    validate_input(K, "K",  must_be_integer = TRUE, not_na = TRUE)
    validate_input(lower_bound, "lower_bound", greater_than = 0, must_be_integer = TRUE)
    validate_input(upper_bound, "upper_bound", greater_than = 0, must_be_integer = TRUE)
    validate_input(beta_max, "beta_max", greater_than = 0, must_be_integer = TRUE)
    validate_input(beta_min, "beta_min", greater_than = -1, must_be_integer = TRUE)
    validate_input(theta_max, "theta_max", greater_than = 0.0)
    validate_input(theta_min, "theta_min", greater_than = 0.0)
    validate_input(phi, "phi", greater_than = 0, smaller_than = 1)
    validate_input(tau, "tau", greater_than = 0, must_be_integer = TRUE)
    validate_input(kmax, "kmax", greater_than = 0, must_be_integer = TRUE)
    validate_input(alpha, "alpha", greater_than = 0.0)
}