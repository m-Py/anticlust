

#' Validating the arguments passed to `anticlustering`
#' 
#' @return NULL
#'
#' @noRd
input_validation_anticlustering <- function(x, K, objective, method,
                                          preclustering, categories,
                                          repetitions, standardize = FALSE, cannot_link = NULL,
                                          must_link = NULL, blocks = NULL) {
  
  ## Validate feature input
  validate_data_matrix(x)
  
  is_df <- inherits(x, "data.frame")
  if (is_df && !is.function(objective)) {
    x <- get_anticlustering_features(x, objective, standardize) # this is the input that must be considered
  }
  
  x <- as.matrix(x)
  N <- nrow(x)
  unequal_group_sizes <- (length(K) != 1) && (sd(table(initialize_clusters(N, K, NULL))) != 0)
  initial_grouping_passed <- length(K) == N

  validate_input(
    method, "method", len = 1,
    input_set = c("ilp", "exchange", "heuristic", "local-maximum", "brusco", "2PML", "3phase"), 
    not_na = TRUE, not_function = TRUE
  )
  
  if (method == "3phase") {
    if (is.function("objective")) stop("objective can only be a function with method = 'exchange' or method = 'local-maximum'.")
    if (objective %in% "dispersion") {
      stop("objective = dispersion does not work with the 3 phase search algorithm.")
    }
    if (argument_exists(categories)) {
      stop("Using the argument categories does not work with the 3 phase search algorithm.")
    }
    if (unequal_group_sizes && objective %in% c("average-diversity", "variance", "kplus")) {
      stop("The three phase algorithm does not support the average diversity, variance and k-plus objective for unequal-sized groups.")
    }
    if (initial_grouping_passed) {
      stop("The 3 phase algorithm does not handle passing an initial grouping via argument `K`.")
    }
  }
  
  if (argument_exists(must_link)) {
    validate_input(must_link, "must_link", not_function = TRUE, len = N)
    must_link <- as.matrix(must_link)

    if (argument_exists(categories)) {
      stop("\nCombining the `categories` argument together with must-link constraints \n",
           "is currently not supported; use the categorical variables as part of the first argument\n",
           "`x` instead (see the vignette on categorical variables).")
    }
    
    validate_input(objective, "objective", input_set = c("diversity","kplus", "variance"), not_na = TRUE, len = 1) 
    validate_input(method, "method", input_set = c("local-maximum", "exchange", "2PML"), not_na = TRUE, len = 1) 
    
    if (unequal_group_sizes && objective %in% c("kplus", "variance")) {
      stop("K-means and k-plus anticlustering are only possible with equal-sized groups when must-link constraints are used.")
    }
    
    if (ncol(must_link) != 1) {
      stop("Argument must_link must be a vector.")
    }

    if (isTRUE(preclustering)) {
      stop("It is not possible to combine preclustering with must-link constraints.")
    }
    if (argument_exists(cannot_link)) {
      stop("Currently, it is not possible to use both cannot-link and must-link constraints.")
    }
  }
  
  if (argument_exists(cannot_link)) {
    cannot_link <- as.matrix(cannot_link)
    validate_input(cannot_link, "cannot_link", 
                   objmode = "numeric", must_be_integer = TRUE, 
                   greater_than = 0, smaller_than = NROW(x)+1, 
                   not_na = TRUE, not_function = TRUE)
    if (ncol(cannot_link) > 2) {
      stop("Argument cannot_link must have 2 columns.")
    }
    if (objective == "dispersion") {
      stop("objective = dispersion does not work with cannot_link constraints.")
    }
    if (!method %in% c("brusco", "exchange", "local-maximum", "ilp")) {
      stop("Cannot-link constraints currently work with method = 'brusco', 'exchange', 'local-maximum', or 'ilp'.")
    }
    if (argument_exists(categories)) {
      stop("\nCombining the `categories` argument together with cannot-link constraints \n",
           "is currently not supported; use the categorical variables as part of the first argument\n",
           "`x` instead (see the vignette on categorical variables).")
    }
    if (isTRUE(preclustering)) {
      stop("It is not possible to combine preclustering with cannot-link constraints.")
    }
  }
  
  validate_input(standardize, "standardize", objmode = "logical", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  
  if (argument_exists(repetitions)) {
    validate_input(repetitions, "repetitions", objmode = "numeric", len = 1, 
                   greater_than = 0, must_be_integer = TRUE, not_na = TRUE,
                   not_function = TRUE)
  }
  
  ## Merge categories variable so that `length` can be applied:
  if (argument_exists(categories)) {
    categories <- merge_into_one_variable(categories)
    validate_input(categories, "categories", not_function = TRUE, len = N)
  }
  if (argument_exists(blocks)) {
    blocks <- merge_into_one_variable(blocks)
    validate_input(blocks, "blocks", not_function = TRUE, len = N)
    if (argument_exists(must_link)) stop("Currently, using must-link constraints is not possible within blocks.")
    if (!method %in% c("exchange", "local-maximum")) {
      stop("Blocked anticlustering currently only works with method = 'exchange' or 'local-maximum'.")
    }
    if (argument_exists(repetitions)) {
      stop("Cannot use multiple repetitions with blocks.")
    }
    if (preclustering) {
      stop("Cannot use preclustering with blocks.")
    }
  }

  validate_input(preclustering, "preclustering", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)

  if (any(is.na(x)) && preclustering) {
    stop("Cannot use preclustering when there is NAs in the input.")
  }
  
  if (method == "2PML") {
    if (!argument_exists(must_link)) {
      stop("Method 2PML only works with must-link constraints.")
    }
  }
  
  if (method == "brusco") {
    if (argument_exists(categories)) {
      stop("It is not possible to use the algorithm by Brusco et al. with categorical restrictions.")
    }
    if (preclustering == TRUE) {
      stop("It is not possible to use the algorithm by Brusco et al. with preclustering restrictions.")
    }
  }
  
                 
  # Allow that K is an initial assignment of elements to clusters
  validate_input(K, "K", objmode = "numeric", must_be_integer = TRUE, not_na = TRUE, not_function = TRUE)
  if (length(K) == 1) {
    validate_input(K, "K", greater_than = 1, smaller_than = N)
  } else {
    if (length(K) != N && sum(K) != N) {
      stop("Argument `K` is misspecified.")
    }
    if (method == "ilp") {
      stop("Pass the number of groups via argument `K` when using method = 'ilp'.")
    }
    if (preclustering == TRUE) {
      stop("Cannot only use preclustering when argument `K` is the number of groups.")
    }
    if (method == "ilp") {
      stop("The argument `K` should indicate the number of groups when using `method = 'ilp'`.")
    }
    if (argument_exists(categories)) {
      if (length(categories) != length(K) && sum(K) != N) {
        stop("Length of arguments `categories` and `K` differ.")
      }
    }
  }

  if (argument_exists(categories) && length(categories) != N) {
    stop("The length of the `categories` argument is not equal to the the number of input elements.")
  }

  if (length(K) == 1 && N %% K != 0) {
    if (method == "ilp") {
      stop("K must be a divider of the number of elements when using the ILP method. ",
           "(Try out method = 'exchange'.)")
    }
  }



  if (method == "ilp") {
    check_if_solver_is_available()
    if (!objective %in% c("distance", "diversity", "dispersion")) {
      stop("The ILP method is only available for the diversity and dispersion objectives.")
    }
    if (argument_exists(repetitions)) {
      stop("Do not use argument `repetitions` when using the ILP method.")
    }
  }

  if (!inherits(objective, "function")) {
    validate_input(objective, "objective", input_set = c(
      "distance", 
      "diversity", 
      "average-diversity", 
      "dispersion",
      "variance", 
      "kplus"
      ), 
      len = 1, not_na = TRUE)
    if (objective %in% c("variance", "kplus") && method == "ilp") {
      stop("You cannot use this function to optimally maximize the kmeans or kplus criterion. Use optimal_anticlustering().")
    }
    if (objective %in% c("variance", "kplus") && is_distance_matrix(x) && !is_df) {
        stop("You cannot use a distance matrix with the objective 'variance' or 'kplus'.")
    }
  }

  if (argument_exists(categories) && method == "ilp") {
    stop("The ILP method cannot incorporate categorical restrictions.")
  }
  return(invisible(NULL))
}


#' A function for input validation
#'
#' @param obj The object that undergoes validation
#' @param argument_name A string indicating the name of the object
#'   This name is used when an error is thrown so the user
#'   is informed on the cause of the error.
#' @param len Optional numeric vector for objects having a length
#'   (mostly for vectors). Tests via `NROW`, so can also test matrix-like 
#'   objects.
#' @param greater_than Optional scalar indicating if numeric input has
#'   to be greater than a specified number.
#' @param must_be_integer Optional logical vector indicating if numeric
#'   input has to be integer.
#' @param groupsize Optional argument how many groups a grouping variable
#'   consist of.
#' @param input_set Optional argument specifying a set of values an
#'   argument can take.
#' @param objmode The required mode of \code{obj}
#' @param not_na Boolean to indicate whether NA input is forbidden
#'   (TRUE means that NA is not allowed)
#' @param smaller_than A value for numeric input that must not be exceeded
#'
#' @return NULL
#'
#' @noRd

validate_input <- function(obj, argument_name, len = NULL, greater_than = NULL,
                           must_be_integer = FALSE, groupsize = NULL, input_set = NULL, 
                           objmode = NULL, not_na = FALSE, smaller_than = NULL, 
                           not_function = NULL) {
                           
                           
  self_validation(argument_name, len, greater_than,
                  must_be_integer, groupsize, input_set,
                  objmode, not_na, not_function)

  argument_name <- paste("Argument", argument_name)
  
  if (argument_exists(not_function) && not_function == TRUE) {
    if (inherits(obj, "function")) {
      stop(argument_name, "must not be a function")
    }
  }
   
  
  ## - Check length of input
  if (argument_exists(len)) {
    if (NROW(obj) != len) {
      stop(argument_name, " must have length ", len)
    }
  }

  ## - Check mode of input
  if (argument_exists(objmode)) {
    if (mode(obj) != objmode) {
      stop(argument_name, " must be ", objmode, ", but is ", mode(obj))
    }
  }
  
  ## - Check if input has to be greater than some value
  if (argument_exists(greater_than)) {
    if (any(obj <= greater_than)) {
      stop(argument_name, " must be greater than ", greater_than)
    }
  }
  
  ## - Check if input has to be greater than some value
  if (argument_exists(smaller_than)) {
    if (any(obj >= smaller_than)) {
      stop(argument_name, " must be smaller than ", smaller_than)
    }
  }
  
  if (not_na == TRUE) {
    if (sum(is.na(obj) >= 1)) {
      stop(argument_name, " must not be NA, but contains NA")
    }
  }
  
  ## - Check if input has to be integer
  if (must_be_integer == TRUE && any(obj %% 1 != 0)) {
    stop(argument_name, " must be integer")
  }
  ## - Check if correct number of groups is provided
  if (argument_exists(groupsize)) {
    if (length(table(obj)[table(obj) != 0]) != groupsize) {
      stop(argument_name, " must consist of exactly ", groupsize,
           " groups with more than 0 observations.")
    }
  }

  ## - Check if argument matches a predefined input set
  if (argument_exists(input_set)) {
    if (!obj %in% input_set) {
      stop(argument_name, " can either be set to '",
           paste(input_set, collapse = "' or '"), "'")
    }
  }

  return(invisible(NULL))
}

## Validate feature input
validate_data_matrix <- function(x, NA_allowed = TRUE) {
  x <- data.frame(x)
  modes <- sapply(x, mode)
  if (any(modes != "numeric")) {
    stop("Your data (the first argument `x`) should only consist of numeric and factor variables.")
  }
  if (any(is.na(x)) && !NA_allowed) {
    stop("Missing values are not allowed.")
  }
}

## Validate input for the `validate_input` function (these errors are
## not for users, but only for developers)
self_validation <- function(argument_name, len, greater_than,
                            must_be_integer, groupsize,
                            input_set, objmode, not_na, not_function) {

  if (argument_exists(not_function)) {
    stopifnot(not_function %in% c(TRUE, FALSE))
  }
                            
  if (argument_exists(len)) {
    stopifnot(class(len) %in% c("numeric", "integer"))
    stopifnot(length(len) == 1)
    stopifnot(len >= 0)
    stopifnot(len %% 1 == 0)
  }

  if (argument_exists(greater_than)) {
    stopifnot(length(greater_than) == 1)
    stopifnot(class(greater_than) %in% c("numeric", "integer"))
  }

  stopifnot(length(must_be_integer) == 1)
  stopifnot(must_be_integer %in% c(TRUE, FALSE))
  stopifnot(length(not_na) == 1)
  stopifnot(not_na %in% c(TRUE, FALSE))

  if (argument_exists(groupsize)) {
    stopifnot(mode(groupsize) == "numeric")
    stopifnot(length(groupsize) == 1)
  }

  if ((argument_exists(input_set) && !argument_exists(len)) ||
      (argument_exists(input_set) && len != 1))  {
    stop("If an input set is passed, argument len must be 1 ",
         "(this message should not be seen by users of the package).")
  }

  if (argument_exists(objmode)) {
    stopifnot(mode(objmode) == "character")
    stopifnot(length(objmode) == 1)
  }

  return(invisible(NULL))
}

argument_exists <- function(arg) {
  !is.null(arg)
}


# Test that the number of features is a multiplier of unique
# anticlusters. I think this function is only used in test cases.
legal_number_of_clusters <- function(features, clusters) {
  ## 1. correct number of clusters assignments?
  if (length(clusters) != nrow(features))
    stop("The number of cluster assignments and the number of features have to match.")
  n_anticlusters <- length(unique(clusters))
  ## 2. More than 1 disctinct cluster?
  if (n_anticlusters <= 1)
    stop("There have to be at least two different anticlusters")
  ## 3. Do all clusters occur equally often?
  if (!all(table(clusters) == table(clusters)[1]))
    stop("Each clusters must occur equally often")
  ## 4. Probably redundant to 3:
  if (nrow(features) %% n_anticlusters != 0)
    stop("The number of elements is not a multiplier of the number of anticlusters")
  invisible(NULL)
}


input_validation_matching <- function(
  x, p, 
  match_between, 
  match_within, 
  match_extreme_first, 
  target_group
) {
  validate_data_matrix(x, FALSE)
  N <- nrow(as.matrix(x))
  validate_input(
    p, "p", 
    len = 1, 
    must_be_integer = TRUE, 
    not_na = TRUE
  )
  if (argument_exists(match_between)) {
    validate_input(
      match_between, 
      "match_between", 
      len = N
    )
  }
  if (argument_exists(match_within)) {
    validate_input(
      match_within, 
      "match_within", 
      len = N
    )
  }
  validate_input(
    match_extreme_first, 
    "match_extreme_first", 
    len = 1, 
    input_set = c(TRUE, FALSE), 
    not_na = TRUE
  )
  if (argument_exists(target_group)) {
    validate_input(
      target_group, 
      "target_group", 
      len = 1, 
      input_set = c("smallest", "diverse", "none"), 
      not_na = TRUE
    )
  }
}

# Check if a solver package can be used
check_if_solver_is_available <- function() {
  glpk_available <- requireNamespace("Rglpk", quietly = TRUE)
  symphony_available <- requireNamespace("Rsymphony", quietly = TRUE)
  lpSolve_available <- requireNamespace("lpSolve", quietly = TRUE)
  gurobi_available <- requireNamespace("gurobi", quietly = TRUE)
  no_solver_available <- !glpk_available && !symphony_available && !lpSolve_available && !gurobi_available
  
  if (no_solver_available) {
    stop("\n\nAn exact method was requested, but no ILP solver is ",
         "available. For example, you could install the GNU linear programming kit: \n\n",
         "- On windows, visit ",
         "http://gnuwin32.sourceforge.net/packages/glpk.htm \n\n",
         "- Use homebrew to install it on mac, 'brew install glpk' \n\n",
         "- 'sudo apt install libglpk-dev' on Ubuntu ",
         "\n\nThen, install the Rglpk package via ",
         "`install.packages('Rglpk')`.\n \n Another possibilty is to install the R package 'Rsymphony' ",
         "and the SYMPHONY ILP solver (https://github.com/coin-or/SYMPHONY), or the R package `lpSolve`.")
  }
  return(invisible(NULL))
}
