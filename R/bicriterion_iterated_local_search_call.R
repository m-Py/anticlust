
#' Bicriterion iterated local search heuristic
#'
#' This function implements the bicriterion algorithm BILS for
#' anticlustering by Brusco et al. (2020;
#' <doi:10.1111/bmsp.12186>). The description of their algorithm is
#' given in Section 3 of their paper (in particular, see the
#' Pseudocode in their Figure 2). As of anticlust version 0.8.6, this
#' function also includes some extensions to the BILS algorithm that
#' are implemented through the optional arguments
#' \code{dispersion_distances}, \code{average_diversity},
#' \code{init_partitions}, and \code{return}. If these arguments are
#' not changed, the function performs the "vanilla" BILS as described
#' in Brusco et al.
#' 
#' 
#' @param x The data input. Can be one of two structures: (1) A
#'     feature matrix where rows correspond to elements and columns
#'     correspond to variables (a single numeric variable can be
#'     passed as a vector). (2) An N x N matrix dissimilarity matrix;
#'     can be an object of class \code{dist} (e.g., returned by
#'     \code{\link{dist}} or \code{\link{as.dist}}) or a \code{matrix}
#'     where the entries of the upper and lower triangular matrix
#'     represent pairwise dissimilarities.
#' @param K How many anticlusters should be created. Alternatively:
#'     (a) A vector describing the size of each group, or (b) a vector
#'     of length \code{nrow(x)} describing how elements are assigned
#'     to anticlusters before the optimization starts.
#' @param R The desired number of restarts for the algorithm. By
#'     default, both phases (MBPI + ILS) of the algorithm are
#'     performed once.  See details.
#' @param W Optional argument, a vector of weights defining the
#'     relative importance of dispersion and diversity (0 <= W <=
#'     1). See details.
#' @param Xi Optional argument, specifies probability of swapping
#'     elements during the iterated local search. See examples.
#' @param dispersion_distances A distance matrix used to compute the
#'     dispersion if the dispersion should not be computed on the
#'     basis of argument \code{x}.
#' @param average_diversity Boolean. Compute the diversity not as a
#'     global sum across all pairwise within-group distances, but as
#'     the sum of the average of within-group distances.
#' @param init_partitions A matrix of initial partitions (rows =
#'     partitions; columns = elements) that serve as starting
#'     partitions during the iterations of the first phase of the BILS
#'     (i.e., the MBPI).  If not passed, a new random partition is
#'     generated at the start of each iteration (which is the default
#'     behaviour).
#' @param return Either "paretoset" (default), "best-diversity", 
#'     "best-average-diversity", "best-dispersion". See below.
#'     
#' @details
#'
#' The bicriterion algorithm by Brusco et al. (2020) aims to
#' simultaneously optimize two anticlustering criteria: the
#' \code{\link{diversity_objective}} and the
#' \code{\link{dispersion_objective}}. It returns a list of partitions
#' that approximate the pareto set of efficient solutions across both
#' criteria. By considering both the diversity and dispersion, this
#' algorithm is well-suited for maximizing overall within-group
#' heterogeneity. To select a partition among the approximated pareto
#' set, it is reasonable to plot the objectives for each partition
#' (see Examples).
#'
#' The arguments \code{R}, \code{W} and \code{Xi} are named for
#' consistency with Brusco et al. (2020). The argument \code{K} is
#' used for consistency with other functions in anticlust; Brusco et
#' al. used `G` to denote the number of groups. However, note that
#' \code{K} can not only be used to denote the number of equal-sized
#' groups, but also to specify group sizes, as in
#' \code{\link{anticlustering}}.
#' 
#' This function implements the combined bicriterion algorithm BILS,
#' which consists of two phases: The multistart bicriterion pairwise
#' interchange heuristic (MBPI, which is a local maximum search
#' heuristic similar to \code{method = "local-maximum"} in
#' \code{\link{anticlustering}}), and the iterated local search (ILS),
#' which is an improvement procedure that tries to overcome local
#' optima.  The argument \code{R} denotes the number of restarts of
#' the two phases of the algorithm. If \code{R} has length 1, half of
#' the repetitions perform the first phase MBPI), the other half
#' perform the ILS.  If \code{R} has length 2, the first entry
#' indicates the number of restarts of MBPI the second entry indicates
#' the number of restarts of ILS.  The argument \code{W} denotes the
#' relative weight given to the diversity and dispersion criterion in
#' a given run of the search heuristic. In each run, the a weight is
#' randomly selected from the vector \code{W}. As default values, we
#' use the weights that Brusco et al. used in their analyses. All
#' values in \code{W} have to be in [0, 1]; larger values indicate
#' that diversity is more important, whereas smaller values indicate
#' that dispersion is more important; \code{w = .5} implies the same
#' weight for both criteria. The argument \code{Xi} is the probability
#' that an element is swapped during the iterated local search
#' (specifically, Xi has to be a vector of length 2, denoting the
#' range of a uniform distribution from which the probability of
#' swapping is selected). For \code{Xi}, the default is selected
#' consistent with the analyses by Brusco et al.
#'
#' If the data input \code{x} is a feature matrix (that is: each row
#' is a "case" and each column is a "variable"), a matrix of the
#' Euclidean distances is computed as input to the algorithm. If a
#' different measure of dissimilarity is preferred, you may pass a
#' self-generated dissimilarity matrix via the argument \code{x}.  The
#' argument \code{dispersion_distances} can additionally be used if
#' the dispersion should be computed on the basis of a different
#' distance matrix.
#' 
#' If multiple \code{init_partitions} are given, ensure that each
#' partition (i.e., each row of\code{init_partitions}) has the exact
#' same output of \code{\link{table}}.
#' 
#' @return By default, a \code{matrix} of anticlustering partitions
#'     (i.e., the approximated pareto set). Each row corresponds to a
#'     partition, each column corresponds to an input element. If the
#'     argument \code{return} is set to either "best-diversity", 
#'     "best-average-diversity", or "best-dispersion", it only returns one 
#'     partition (as a vector), that maximizes the respective objective.
#' 
#' @author Martin Breuer \email{M.Breuer@@hhu.de}, Martin Papenberg
#'     \email{martin.papenberg@@hhu.de}
#' 
#' @export
#' 
#' @examples 
#' # Generate some random data
#' M <- 3
#' N <- 80
#' K <- 4
#' data <- matrix(rnorm(N * M), ncol = M)
#'
#' # Perform bicriterion algorithm, use 200 repetitions
#' pareto_set <- bicriterion_anticlustering(data, K = K, R = 200)
#'
#' # Compute objectives for all solutions
#' diversities_pareto <- apply(pareto_set, 1, diversity_objective, x = data)
#' dispersions_pareto <- apply(pareto_set, 1, dispersion_objective, x = data)
#'
#' # Plot the pareto set
#' plot(
#'   diversities_pareto,
#'   dispersions_pareto,
#'   col = "blue",
#'   cex = 2,
#'   pch = as.character(1:NROW(pareto_set))
#' )
#' 
#' # Get some random solutions for comparison
#' rnd_solutions <- t(replicate(n = 200, sample(pareto_set[1, ])))
#' 
#' # Compute objectives for all random solutions
#' diversities_rnd <- apply(rnd_solutions, 1, diversity_objective, x = data)
#' dispersions_rnd <- apply(rnd_solutions, 1, dispersion_objective, x = data)
#' 
#' # Plot random solutions and pareto set. Random solutions are far away 
#' # from the good solutions in pareto set
#' plot(
#'   diversities_rnd, dispersions_rnd, 
#'   col = "red",
#'   xlab = "Diversity",
#'   ylab = "Dispersion",
#'   ylim = c(
#'     min(dispersions_rnd, dispersions_pareto), 
#'     max(dispersions_rnd, dispersions_pareto)
#'   ),
#'   xlim = c(
#'     min(diversities_rnd, diversities_pareto), 
#'     max(diversities_rnd, diversities_pareto)
#'   )
#' )
#' 
#' # Add approximated pareto set from bicriterion algorithm:
#' points(diversities_pareto, dispersions_pareto, col = "blue", cex = 2, pch = 19)
#'
#' @note
#'
#' For technical reasons, the pareto set returned by this function has
#' a limit of 500 partitions. Usually however, the
#' algorithm usually finds much fewer partitions. There is one following exception:
#' We do not recommend to use this method when the input data is
#' one-dimensional where the algorithm may identify too many
#' equivalent partitions causing it to run very slowly (see section 5.6 in 
#' Breuer, 2020).
#' 
#' @references
#' 
#' Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining
#' diversity and dispersion criteria for anticlustering: A bicriterion
#' approach. British Journal of Mathematical and Statistical
#' Psychology, 73, 275-396. https://doi.org/10.1111/bmsp.12186
#' 
#' Breuer (2020). Using anticlustering to maximize diversity and dispersion:
#' Comparing exact and heuristic approaches. Bachelor thesis.
#' 

bicriterion_anticlustering <- function(
  x, K, R = NULL, 
  W = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 0.99, 0.999, 0.999999),
  Xi = c(0.05, 0.1),
  dispersion_distances = NULL, average_diversity = FALSE, init_partitions = NULL, return = "paretoset") {
  
  input_validation_bicriterion_anticlustering(x, K, R, W, Xi, dispersion_distances, average_diversity, init_partitions, return)

  distances <- convert_to_distances(x) 
  N <- NROW(distances)
  WL <- length(W)
  if (is.null(R)) {
    R <- c(1, 1)
  } else if (length(R) == 1) {
    R <- rep(ceiling(R / 2), 2)
  }
  
  if (is.null(dispersion_distances)) {
    dispersion_distances <- distances
  } else {
    dispersion_distances <- convert_to_distances(dispersion_distances)
  }

  clusters <- initialize_clusters(N, K, NULL) - 1
  
  if (argument_exists(init_partitions)) {
    use_init_partitions <- 1
    init_partitions <- t(apply(init_partitions, 1, to_numeric)) - 1 # ensure correct labeling of partitions
    clusters <- init_partitions[1, ]
  } else {
    init_partitions <- 0
    use_init_partitions <- 0
  }
  
  frequencies <- rep(1, length(unique(clusters)))
  if (average_diversity) {
    frequencies <- table(clusters)
  }
  upper_bound = 500 # limits number of resulting partitions 
  # create empty matrix for results to use in C
  result_matrix = matrix(data = -1, nrow = upper_bound , ncol = N) 
  
  # Call C function
  results <- .C(
    "bicriterion_iterated_local_search_call",
    as.double(distances),
    as.double(dispersion_distances),
    as.integer(N),
    as.integer(R),
    as.integer(upper_bound),
    as.integer(WL),
    as.double(W),
    as.double(Xi),
    as.integer(clusters),
    as.integer(frequencies),
    as.integer(use_init_partitions),
    as.integer(t(init_partitions)),
    result = integer(length(result_matrix)),
    mem_error = as.integer(0),
    PACKAGE = "anticlust" # important to call C
  )

  if (results[["mem_error"]] == 1) {
    stop("Could not allocate enough memory.")
  }
  
  results <- results$result
  
  # remove allocated space that was not needed
  results <- results[results != -1]

  # C returns the list of partitions as one vector, 
  # that we turn back into a matrix
  results <- matrix(results, ncol = N, byrow = TRUE)
  # Order each partition:
  results <- data.frame(t(apply(results, 1, order_cluster_vector)))
  # Remove duplicates
  results <- results[!duplicated(results), ]
  results <- unname(as.matrix(results))
  if (return == "paretoset") {
    return(results)
  }
  if (return == "best-dispersion") {
    best_obj <- which.max(apply(results, 1, dispersion_objective_, dispersion_distances))
    return(results[best_obj, ])
  } else if (return == "best-average-diversity" || average_diversity) {
    # manually compute objectives, apply() does not work here because table() has different results for the different clusterings
    average_diversities <- c()
    for (i in 1:nrow(results)) {
      average_diversities[i] <- weighted_diversity_objective_(distances, results[i, ], table(results[i, ]))
    }
    return(results[which.max(average_diversities), ])
  } else if (return == "best-diversity") {
    best_obj <- which.max(apply(results, 1, diversity_objective_, distances))
    return(results[best_obj, ])
  } 
}

input_validation_bicriterion_anticlustering <- function(
    x, K, R, W, Xi, dispersion_distances, 
    average_diversity, init_partitions, return) {

  input_validation_anticlustering(
    x, K, objective = "diversity", method = "brusco", 
    preclustering = FALSE, categories = NULL,
    repetitions = 1, standardize = FALSE
  )
  
  x <- convert_to_distances(x)
  N <- nrow(x)
  
  validate_input(
    return, "return", objmode = "character", len = 1,
    input_set = c("paretoset", "best-diversity", "best-average-diversity", "best-dispersion"), 
    not_na = TRUE, 
    not_function = TRUE
  )
  
  validate_input(
    average_diversity, "average_diversity", objmode = "logical", len = 1,
    input_set = c(TRUE, FALSE), 
    not_na = TRUE, 
    not_function = TRUE
  )
  
  if (average_diversity && return == "best-diversity") {
    stop("Cannot use 'return = best-diversity' with 'average_diversity = TRUE' (probably you want to use 'return = best-average-diversity').")
  }
  if (!average_diversity && return == "best-average-diversity") {
    stop("Cannot use 'return = best-average-diversity' with 'average_diversity = FALSE' (probably you want to use 'return = best-diversity').")
  }

  if (argument_exists(dispersion_distances)) {
    if (any(dim(dispersion_distances) != dim(x))) {
      stop("The dimensions of argument 'x' and 'dispersion_distances' do not match.")
    } 
  }

  checkweights(W)
  checkneighborhood(Xi)
  
  if (argument_exists(R)) {
    validate_input(R, "R", must_be_integer = TRUE, not_na = TRUE, not_function = TRUE, greater_than = -1)
    if (!length(R) %in% 1:2) {
      stop("Argument 'R' must have length 1 or 2.")
    }
  }
  
  if (argument_exists(init_partitions)) {
    if (length(R) != 2) {
      stop("If 'init_partitions' is used, specify argument R as vector of length 2.")
    }
    if (nrow(init_partitions) != R[1]) {
      stop("The number of repetitions (argument 'R') and the number of partitions (argument 'init_partitions') do not match.")
    } 
    if (ncol(init_partitions) != N) {
      stop("The number of columns in 'init_partitions' does not match the number of cases in 'x'.")
    }
  }
}


#verify input
#############
checkweights <- function(W){
  validate_input(W, "W", not_na = TRUE, not_function = TRUE, objmode = "numeric")
  for(i in W){
    if(i < 0 | i > 1){
      stop("All weights (argument 'W') must be percentages. Only values between 0 and 1 are allowed.")
    }
  }  
}

checkneighborhood <- function(Xi) {
  validate_input(Xi, "Xi", not_na = TRUE, not_function = TRUE, len = 2, objmode = "numeric")
  for (i in Xi) {
    if (i < 0 | i > 1) {
      stop("Both neighorhood indexes must be a percentage. Only values between 0 and 1 are allowed.")
    }
  } 
  if (Xi[1] > Xi[2]) {
    stop("First neighborhood percentage needs to be smaller than the second.")
  }
}


