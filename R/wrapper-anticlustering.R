
#' Anticlustering
#'
#' Partition a pool of elements into groups (i.e., anticlusters) with
#' the aim of creating high within-group heterogeneity and high
#' between-group similarity.  Anticlustering is accomplished by
#' maximizing instead of minimizing a clustering objective function.
#' Implements anticlustering methods as described in Papenberg and
#' Klau (2021; <doi:10.1037/met0000301>), Brusco et al. 
#' (2020; <doi:10.1111/bmsp.12186>), Papenberg (2024; 
#' <doi:10.1111/bmsp.12315>), and Papenberg et al. (2025; 
#' <doi:10.1101/2025.03.03.641320>).
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
#' @param objective The objective to be maximized. The options
#'     "diversity" (default; previously called "distance", which is
#'     still supported), "average-diversity", "variance", "kplus" and "dispersion" are
#'     natively supported. May also be a user-defined function. See
#'     Details.
#' @param method One of "exchange" (default) , "local-maximum",
#'     "brusco", "ilp", or "2PML".  See Details.
#' @param preclustering Boolean. Should a preclustering be conducted
#'     before anticlusters are created? Defaults to \code{FALSE}. See
#'     Details.
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical variables whose distribution should be similar 
#'     between groups. See Details.
#' @param repetitions The number of times a search heuristic is
#'     initiated when using \code{method = "exchange"}, \code{method =
#'     "local-maximum"}, \code{method = "brusco"}, or \code{method = "2PML"}. 
#'     In the end, the best objective found across the repetitions is returned.
#' @param standardize Boolean. If \code{TRUE} and \code{x} is a
#'     feature matrix, the data is standardized through a call to
#'     \code{\link{scale}} before the optimization starts. This
#'     argument is silently ignored if \code{x} is a distance matrix.
#' @param cannot_link A 2 column matrix where each row has the indices 
#'     of two elements that must not be assigned to the same anticluster. 
#'     Alternatively a vector of length \code{nrow(x)}; elements having 
#'     the same value in this vector cannot be assigned to the same anticluster.
#' @param must_link A numeric vector of length \code{nrow(x)}. Elements having 
#'     the same value in this vector are assigned to the same anticluster.
#'
#' @return A vector of length N that assigns a group (i.e, a number
#'     between 1 and \code{K}) to each input element.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist dist sd
#' 
#' @useDynLib anticlust, .registration = TRUE
#'
#' @export
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @details
#'
#' This function is used to solve anticlustering. That is, the data
#' input is divided into \code{K} groups in such a way that elements
#' within groups are heterogeneous and the different groups are
#' similar. Anticlustering is accomplished by maximizing instead of
#' minimizing a clustering objective function. The maximization of
#' five objectives is natively supported (other
#' functions can also defined by the user as described below):
#' 
#' \itemize{
#'   \item{the diversity, setting 
#'         \code{objective = "diversity"} (this is the default objective)}
#'   \item{the average diversity, which normalizes the diversity by cluster size, 
#'        setting \code{objective = "average-diversity"}}
#'   \item{the k-means (or "variance") objective, setting \code{objective = "variance"}}
#'   \item{the k-plus objective, an extension of the k-means objective,
#'         setting \code{objective = "kplus"}}
#'   \item{the dispersion, which is the minimum distance between 
#'         any two elements within the same cluster (setting 
#'         \code{objective = "dispersion"})}
#' }
#'
#' The k-means objective is the within-group variance---that is, the
#' sum of the squared distances between each element and its cluster
#' center (see \code{\link{variance_objective}}). K-means
#' anticlustering focuses on minimizing differences with regard to the
#' means of the input variables (that is, the columns in \code{x}), but it ignores any other distribution
#' characterstics such as the variance / standard deviation. K-plus anticlustering
#' (using \code{objective = "kplus"}) is an extension of the k-means criterion that also
#' minimizes differences with regard to the standard
#' deviations between groups (for details see \code{\link{kplus_anticlustering}}). K-plus
#' anticlustering can also be extended towards higher order moments such as skew and kurtosis; 
#' to consider these additional distribution characteristics, use the function
#' \code{\link{kplus_anticlustering}}. Setting \code{objective = "kplus"} in 
#' \code{anticlustering} function will only consider means 
#' and standard deviations (in my experience, this is what users usually want). 
#' It is strongly recommended to set the argument \code{standardize = TRUE} 
#' when using the k-plus objective.
#' 
#' The "diversity" objective is the sum of pairwise
#' distances of elements within the same groups (see
#' \code{\link{diversity_objective}}). Hence, anticlustering using the diversity 
#' criterion maximizes between-group similarity
#' by maximizing within-group heterogeneity (represented as the sum of all pairwise distances). 
#' If it is computed on the basis of the Euclidean distance (which is the default
#' behaviour if \code{x} is a feature matrix), the diversity is an all rounder objective that 
#' tends to equate all distribution 
#' characteristics between groups (such as means, variances, ...). 
#' Note that the equivalence of within-group heterogeneity and between-group similarity only
#' holds for equal-sized groups. For unequal-sized groups, it is recommended to
#' use a different objective when striving for overall between-group similarity,
#' e.g., the k-plus objective or the \code{"average-diversity"}. The average diversity
#' was introduced in version 0.8.6, and it is more useful if groups are not 
#' equal-sized. The average diversity normalizes the sum of intra-cluster distances 
#' by group size. If all groups are equal-sized, it is equivalent to the 
#' regular diversity. In the publication that introduces
#' the \code{anticlust} package (Papenberg & Klau, 2021), we used the term "anticluster 
#' editing" to refer to the maximization of the diversity, because the reversed 
#' procedure - minimizing the diversity - is also known as "cluster editing". 
#' 
#' The "dispersion" is the minimum distance between any two elements
#' that are part of the same cluster; maximization of this objective
#' ensures that any two elements within the same group are as
#' dissimilar from each other as possible. Applications that require
#' high within-group heterogeneity often require to maximize the
#' dispersion. Oftentimes, it is useful to also consider the diversity
#' and not only the dispersion; to optimize both objectives at the
#' same time, see the function
#' \code{\link{bicriterion_anticlustering}}.
#'
#' If the data input \code{x} is a feature matrix (that is: each row
#' is a "case" and each column is a "variable") and the option
#' \code{objective = "diversity"} or \code{objective = "dispersion"} is used, 
#' the Euclidean distance is computed as the basic unit of the objectives. If
#' a different measure of dissimilarity is preferred, you may pass a
#' self-generated dissimilarity matrix via the argument \code{x}.
#'
#' In the standard case, groups of equal size are generated. Adjust
#' the argument \code{K} to create groups of different size (see
#' Examples).
#'
#' \strong{Algorithms for anticlustering}
#'
#' By default, a heuristic method is employed for anticlustering: the
#' exchange method (\code{method = "exchange"}). First, elements are
#' randomly assigned to anticlusters (It is also possible to
#' explicitly specify the initial assignment using the argument
#' \code{K}; in this case, \code{K} has length \code{nrow(x)}.) Based
#' on the initial assignment, elements are systematically swapped
#' between anticlusters in such a way that each swap improves the
#' objective value. For an element, each possible swap with elements
#' in other clusters is simulated; then, the one swap is performed
#' that improves the objective the most, but a swap is only conducted
#' if there is an improvement at all. This swapping procedure is
#' repeated for each element. When using \code{method =
#' "local-maximum"}, the exchange method does not terminate after the
#' first iteration over all elements; instead, the swapping continues
#' until a local maximum is reached. This method corresponds to the algorithm 
#' "LCW" by Weitz and Lakshminarayanan (1998). This means that after the
#' exchange process has been conducted once for each data point, the
#' algorithm restarts with the first element and proceeds to conduct
#' exchanges until the objective cannot be improved.
#'
#' When setting \code{preclustering = TRUE}, only the \code{K - 1}
#' most similar elements serve as exchange partners for each element,
#' which can speed up the optimization (more information
#' on the preclustering heuristic follows below). If the \code{categories} argument
#' is used, only elements having the same value in \code{categories} serve as exchange
#' partners.
#' 
#' Using \code{method = "brusco"} implements the local bicriterion
#' iterated local search (BILS) heuristic by Brusco et al. (2020) and
#' returns the partition that best optimized either the diversity or
#' the dispersion during the optimization process. The function
#' \code{\link{bicriterion_anticlustering}} can also be used to run
#' the algorithm by Brusco et al., but it returns multiple partitions
#' that approximate the optimal pareto efficient set according to both
#' objectives (diversity and dispersion). Thus, to fully utilize the
#' BILS algorithm, use the function
#' \code{\link{bicriterion_anticlustering}}.
#'
#' \strong{Optimal anticlustering}
#'
#' Usually, heuristics are employed to tackle anticlustering problems,
#' and their performance is generally very satisfying.  However,
#' heuristics do not investigate all possible group assignments and
#' therefore do not (necessarily) find the
#' "globally optimal solution", i.e., a partitioning that has the best
#' possible value with regard to the objective that is optimized.  Enumerating
#' all possible partitions in order to find the best solution,
#' however, quickly becomes impossible with increasing N, and
#' therefore it is not possible to find a global optimum this
#' way. Because all anticlustering problems considered here are also
#' NP-hard, there is also no (known) clever algorithm that might
#' identify the best solution without considering all possibilities -
#' at least in the worst case. Integer linear programming (ILP) is an
#' approach for tackling NP hard problems that nevertheless tries to
#' be clever when finding optimal solutions: It does not necessarily
#' enumerate all possibilities but is still guaranteed to return the
#' optimal solution. Still, for NP hard problems such as
#' anticlustering, ILP methods will also fail at some point (i.e.,
#' when N increases).
#'
#' \code{anticlust} implements optimal solution algorithms via integer
#' linear programming. In order to use the ILP methods, set
#' \code{method = "ilp"}. The integer linear program optimizing the
#' diversity was described in Papenberg & Klau, (2021; (8) -
#' (12)). It can also be used to optimize the k-means and k-plus objectives,
#' but you actually have to use the function \code{\link{optimal_anticlustering}}
#' for these objectives. The documentation of the function
#' \code{\link{optimal_dispersion}} and \code{\link{optimal_anticlustering}}
#' contain more information on the optimal anticlustering algorithms.
#' 
#' \strong{Categorical variables}
#'
#' There are two ways to balance categorical variables among anticlusters (also 
#' see the package vignette "Using categorical variables with anticlustering").
#' The first way is to treat them as "hard constraints" via the argument 
#' \code{categories} (see Papenberg & Klau, 2021). If done so, balancing the 
#' categorical variable is accomplished via \code{categorical_sampling} through
#' a stratified split before the anticlustering optimization. After that, the 
#' balance is never changed when the algorithm runs (hence, it is a "hard constraint"). 
#' When \code{categories} has multiple columns (i.e., there are multiple 
#' categorical variables), each combination of categories is treated as a
#'  distinct category by the exchange method (i.e., the multiple columns
#' are "merged" into a single column). This behaviour may lead
#' to less than optimal results on the level of each single categorical variable.
#' In this case, it may be useful to treat the categorical variables as part of 
#' the numeric data, i.e., the first argument \code{x} via binary coding 
#' (e.g. using \code{\link{categories_to_binary}}). The examples show how to do this 
#' when using the bicriterion algorithm by Bruso et al. Using the argument 
#' \code{categories} is only available for the classical exchange procedures, 
#' that is, for \code{method = "exchange"} and \code{method = "local-maximum"}. 
#' 
#' \strong{Anticlustering with constraints}
#' 
#' Versions 0.8.6 and 0.8.7 of anticlust introduced the possibility to induce
#' cannot-link and must-link constraints with anticlustering with 
#' the arguments \code{cannot_link} and \code{must_link}, respectively.
#' Cannot-link constraints ensure that pairs of items are assigned to different
#' clusters. They are given as a 2-column matrix, where each row has the indices
#' of two elements, which must not be assigned to the same cluster. It is possible
#' that a set of cannot-link constraints cannot be fulfilled. To verify whether
#' the constraints cannot be fulfilled (and to actually assign elements 
#' while respecting the constraints), a graph coloring algorithm algorithm is used.
#' This algorithm is is actually the same method as used in \code{\link{optimal_dispersion}}. 
#' The graph coloring algorithm uses an ILP solver and it greatly profits (that is,
#' it may be much faster) from the Rsymphony package, which is not installed as 
#' a necessary dependency with anticlust. It is therefore recommended to 
#' manually install the Rsymphony package, which is then automatically 
#' selected as solver when using the \code{must_link} argument. If you have
#' access to the gurobi solver and have the gurobi R package installed, it will
#' be selected as solver (which is even faster than Symphony). As of version 0.8.11, 
#' it is also possible to use \code{cannot_link} as a vector. In this case it 
#' is ensured that elements having the same value in \code{cannot_link} are not
#' linked in the same cluster. 
#' 
#' Must-link constraints are passed as a single vector of length \code{nrow(x)}.
#' Positions that have the same numeric index are assigned to the same anticluster 
#' (if the constraints can be fulfilled). When including must-link constraints, 
#' \code{method = "2PML"} performs a specialized search heuristic that potentially
#' yields better results than \code{method = "local-maximum"}. The must-link 
#' functionality and the 2PML algorithm was introduced in Papenberg et al. (2025).
#' 
#' The examples illustrate the usage of the \code{must_link} and \code{cannot_link}
#' arguments. Currently, the different kinds of constraints (arguments \code{must_link}, 
#' \code{cannot_link}, and \code{categories}) cannot be used together, but this 
#' may change in future versions. 
#' 
#' \strong{Preclustering}
#' 
#' A useful heuristic for anticlustering is to form small groups of
#' very similar elements and assign these to different groups. This
#' logic is used as a preprocessing when setting \code{preclustering =
#' TRUE}. That is, before the anticlustering objective is optimized, a
#' cluster analysis identifies small groups of similar elements (pairs
#' if \code{K = 2}, triplets if \code{K = 3}, and so forth). The
#' optimization of the anticlustering objective is then conducted
#' under the constraint that these matched elements cannot be assigned
#' to the same group. When using the exchange algorithm, preclustering
#' is conducted using a call to \code{\link{matching}}. When using
#' \code{method = "ilp"}, the preclustering optimally finds groups of
#' minimum pairwise distance by solving the integer linear program
#' described in Papenberg and Klau (2021; (8) - (10), (12) - (13)).
#' Note that when combining preclustering restrictions with \code{method = "ilp"},
#' the anticlustering result is no longer guaranteed to be globally optimal, but
#' only optimal given the preclustering restrictions.
#' 
#' 
#' \strong{Optimize a custom objective function}
#' 
#' It is possible to pass a \code{function} to the argument
#' \code{objective}. See \code{\link{dispersion_objective}} for an
#' example. If \code{objective} is a function, the exchange method
#' assigns elements to anticlusters in such a way that the return
#' value of the custom function is maximized (hence, the function
#' should return larger values when the between-group similarity is
#' higher). The custom function has to take two arguments: the first
#' is the data argument, the second is the clustering assignment. That
#' is, the argument \code{x} will be passed down to the user-defined
#' function as first argument. \strong{However, only after}
#' \code{\link{as.matrix}} has been called on \code{x}. This implies
#' that in the function body, columns of the data set cannot be
#' accessed using \code{data.frame} operations such as
#' \code{$}. Objects of class \code{dist} will be converted to matrix
#' as well.
#' 
#' 
#' @examples
#'
#' # Use default method ("exchange") and the default diversity criterion, also include
#' # a categorical variable via argument `categories`:
#' anticlusters <- anticlustering(
#'   schaper2019[, 3:6],
#'   K = 3,
#'   categories = schaper2019$room
#' )
#' # Compare feature means and standard deviations by anticluster
#' mean_sd_tab(schaper2019[, 3:6], anticlusters)
#' # Verify that the "room" is balanced across anticlusters:
#' table(anticlusters, schaper2019$room)
#' 
#' # Use multiple starts of the algorithm to improve the objective and
#' # optimize the k-means criterion ("variance")
#' anticlusters <- anticlustering(
#'   schaper2019[, 3:6],
#'   objective = "variance",
#'   K = 3,
#'   categories = schaper2019$room,
#'   method = "local-maximum", # better search algorithm
#'   repetitions = 20 # multiple restarts of the algorithm 
#' )
#' # Compare means and standard deviations by anticluster
#' mean_sd_tab(schaper2019[, 3:6], anticlusters)
#' 
#' # Use different group sizes and optimize the extended k-means
#' # criterion ("kplus")
#' anticlusters <- anticlustering(
#'   schaper2019[, 3:6],
#'   objective = "kplus",
#'   K = c(24, 24, 48),
#'   categories = schaper2019$room,
#'   repetitions = 20,
#'   method = "local-maximum",
#'   standardize = TRUE # ususally recommended
#' )
#' 
#' # Use cannot_link constraints: Element 1 must not be linked with elements 2 to 10:
#' cl_matrix <- matrix(c(rep(1, 9), 2:10), ncol = 2)
#' cl <- anticlustering(
#'   schaper2019[, 3:6],
#'   K = 10,
#'   cannot_link = cl_matrix
#' )
#' all(cl[1] != cl[2:10])
#' 
#' # Use cannot_link constraints: Element 1 must be linked with elements 2 to 10.
#' # Element 11 must be linked with elements 12-20.
#' must_link <- rep(NA, nrow(schaper2019))
#' must_link[1:10] <- 1
#' must_link[11:20] <- 2
#' cl <- anticlustering(
#'   schaper2019[, 3:6],
#'   K = 3,
#'   must_link = must_link
#' )
#' cl[1:10]
#' cl[11:20]
#' 
#' # Use the heuristic by Brusco et al. (2020) for k-plus anticlustering
#' # Include categorical variable as part of the optimization criterion rather 
#' # than the argument categories!
#' anticlusters <- anticlustering(
#'   cbind(
#'     kplus_moment_variables(schaper2019[, 3:6], 2), 
#'     categories_to_binary(schaper2019$room)
#'   ),
#'   objective = "variance", # k-plus anticlustering because of the input above!
#'   K = 3,
#'   repetitions = 20,
#'   method = "brusco"
#' )
#' 
#' mean_sd_tab(schaper2019[, 3:6], anticlusters)
#' table(anticlusters, schaper2019$room)
#' 
#' @references
#' 
#' Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining
#' diversity and dispersion criteria for anticlustering: A bicriterion
#' approach. British Journal of Mathematical and Statistical
#' Psychology, 73, 275-396. https://doi.org/10.1111/bmsp.12186
#' 
#' Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition 
#' data sets into equivalent parts. Psychological Methods, 26(2), 
#' 161–174. https://doi.org/10.1037/met0000301.
#' 
#' Papenberg, M. (2024). K-plus Anticlustering: An Improved k-means Criterion for 
#' Maximizing Between-Group Similarity. British Journal of Mathematical and 
#' Statistical Psychology, 77(1), 80-102. https://doi.org/10.1111/bmsp.12315
#' 
#' Papenberg, M., Wang, C., Diop, M., Bukhari, S. H., Oskotsky, B., Davidson, 
#' B. R., Vo, K. C., Liu, B., Irwin, J. C., Combes, A., Gaudilliere, B., 
#' Li, J., Stevenson, D. K., Klau, G. W., Giudice, L. C., Sirota, M., 
#' & Oskotsky, T. T. (2025). Anticlustering for sample allocation to minimize 
#' batch effects. bioRxiv. https://doi.org/10.1101/2025.03.03.641320
#'
#' Späth, H. (1986). Anticlustering: Maximizing the variance criterion.
#' Control and Cybernetics, 15, 213-218.
#' 
#' Weitz, R. R., & Lakshminarayanan, S. (1998). An empirical comparison of 
#' heuristic methods for creating maximally diverse groups. Journal of the 
#' Operational Research Society, 49(6), 635-646. https://doi.org/10.1057/palgrave.jors.2600510
#'

anticlustering <- function(x, K, objective = "diversity", method = "exchange",
                           preclustering = FALSE, categories = NULL, 
                           repetitions = NULL, standardize = FALSE, cannot_link = NULL,
                           must_link = NULL) {


  ## Get data into required format
  input_validation_anticlustering(x, K, objective, method, preclustering, 
                                  categories, repetitions, standardize, cannot_link,
                                  must_link)

  x <- to_matrix(x)
  N <- nrow(x)
  # there is a reason why scaling happens here and below (because of ILP + kplus)
  if (!is_distance_matrix(x) && standardize == TRUE) {
    x <- scale(x)
  }
  
  NUMBER_OF_ANTICLUSTERS <- length(table(initialize_clusters(N, K, NULL)))
  TARGET_GROUPS <- table(initialize_clusters(N, K, NULL))

  ## Exact method using ILP
  if (method == "ilp") {
    if (objective == "dispersion") {
      return(optimal_dispersion(x, K)$groups)
    }
    return(exact_anticlustering(x, K, preclustering, cannot_link))
  }
  
  if (argument_exists(cannot_link)) {
    cannot_link <- as.matrix(cannot_link)
    init <- initialize_clusters(N, K, NULL) # only used if clustering is already given
    if (length(K) != N) {
      if (NCOL(cannot_link) == 1) {# cannot_link is an ID/grouping vector
        if (max(table(c(cannot_link))) > K) {
          stop("Cannot-link constraints cannot be fulfilled.")
        }
        init <- t(replicate(max(1, repetitions), categorical_sampling(c(cannot_link), K = K)))
      } else if (NCOL(cannot_link) == 2) {
        init <- optimal_cannot_link(N, NUMBER_OF_ANTICLUSTERS, table(init), cannot_link, repetitions)
      }
    }

    return(cannot_link_anticlustering(
      x = x, 
      init_clusters = init,
      cannot_link = cannot_link,
      objective = objective,
      method = method
    ))
  }
  
  if (argument_exists(must_link)) {
    return(
      must_link_anticlustering(
        convert_to_distances(x), 
        K, must_link = must_link, 
        method = method, 
        objective = "diversity", 
        repetitions = repetitions
      )
    )
  }

  # Preclustering and categorical constraints are both processed in the
  # variable `categories` after this step:
  categories <- get_categorical_constraints(x, K, preclustering, categories)

  if (!inherits(objective, "function")) {
    if (objective == "kplus") {
      x <- cbind(x, squared_from_mean(x))
      objective <- "variance"
    }
    if (objective == "distance") {
      objective <- "diversity"
    }
  }
  
  if (!is_distance_matrix(x) && standardize == TRUE) {
    x <- scale(x)
  }
  
  # BILS by Brusco et al.:
  if (method == "brusco") {
    average_diversity <- ifelse(objective == "average-diversity", TRUE, FALSE)
    if (objective == "kplus") {
      x <- kplus_moment_variables(x, 2)
      objective <- "variance"
    } 
    if (objective == "variance") {
      x <- convert_to_distances(x)^2
      average_diversity <- TRUE
      objective <- "average-diversity"
    }
    return(bicriterion_anticlustering(x, K, repetitions, average_diversity = average_diversity, return = paste0("best-", objective)))
  }
  
  # Some special cases must be considered now:
  # (a) Is a user defined objective function passed? 
  # (b) Do we need "repeated" anticlustering (i.e., calling the standard exchange method multiple times via
  #     local maximum search and/or multiple initializations)
  # (c) Is the repeated optimization implemented in C, or do we just call the exchange method repeatedly from R?
  
  is_objective_user_defined <- inherits(objective, "function")
  repeated_anticlustering_needed <- method == "local-maximum" || (method == "exchange" && argument_exists(repetitions))
  repeated_anticlustering_has_c_implementation <- !inherits(objective, "function") && objective %in% c("diversity", "average-diversity")
  
  # Diversity anticlustering has C implementation for repeated anticlustering, consider this special case:
  if (repeated_anticlustering_needed && !repeated_anticlustering_has_c_implementation) { 
    repetitions <- ifelse(!argument_exists(repetitions), 1, repetitions)
      return(repeat_anticlustering(x, K, objective, categories, method, repetitions))
  }
  # Most generic exchange method for user defined objectives:
  if (is_objective_user_defined) {
    return(exchange_method(x, K, objective, categories))
  }

  # Redirect to specialized fast exchange methods for diversity, dispersion, kmeans/kplus objectives:
  local_maximum <- ifelse(method == "local-maximum", TRUE, FALSE)
  if (argument_exists(repetitions) && repetitions > 1) { # this can only be the case for diversity objective, based on the logic above
    repetitions <- t(simplify2array(get_multiple_initial_clusters(N, K, categories, repetitions))) - 1
  } else if (argument_exists(repetitions) && repetitions == 1) {
    repetitions <- NULL
  }
  c_anticlustering(x, K, categories, objective, local_maximum = local_maximum, init_partitions = repetitions)
}

# Function that processes input and returns the data set that the
# optimization is conducted on as matrix (for exchange method)
# Returned matrix either represents distances or features.
to_matrix <- function(data) {
  if (!is_distance_matrix(data)) {
    data <- as.matrix(data)
    return(data)
  }
  as.matrix(as.dist(data))
}

# Determines if preclustering constraints or categorical constraints
# are present. Returns a grouping vector if one or both constraints 
# have been passed, or NULL if none is required
get_categorical_constraints <- function(data, K, preclustering, categories) {
  if (preclustering == TRUE) {
    N <- nrow(data)
    matches <- matching(data, p = K, match_within = categories, sort_output = FALSE)
    # deal with NA in matches
    return(replace_na_by_index(matches))
  }
  if (argument_exists(categories)) {
    return(merge_into_one_variable(categories))
  }
  NULL
}

replace_na_by_index <- function(matches) {
  na_matches <- is.na(matches)
  NAs <- sum(na_matches) 
  if (NAs == 0) {
    return(matches)
  }
  max_group <- max(matches, na.rm = TRUE)
  matches[na_matches] <- max_group + 1:NAs 
  matches
}
