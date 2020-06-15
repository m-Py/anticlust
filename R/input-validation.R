

#' Validating the arguments passed to `anticlustering`
#' 
#' @return NULL
#'
#' @noRd
input_validation_anticlustering <- function(x, K, objective, method,
                                          preclustering, categories) {

  ## Merge categories variable so that `length` can be applied:
  categories <- merge_into_one_variable(categories)
  
  
  ## Validate feature input
  validate_data_matrix(x)
  x <- as.matrix(x)
  N <- nrow(x)

  validate_input(preclustering, "preclustering", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE)

  # Allow that K is an initial assignment of elements to clusters
  if (length(K) == 1) {
    validate_input(K, "K", objmode = "numeric", len = 1, 
                   greater_than = 1, must_be_integer = TRUE, not_na = TRUE)
  } else {
    validate_input(K, "K", objmode = "numeric", len = N, not_na = TRUE)
    if (method != "exchange") {
      stop("Passing an initial cluster assignment via the argument `K` ",
           "only works with method = 'exchange'")
    }
    if (preclustering == TRUE) {
      stop("Cannot use preclustering because elements have already been assigned ",
           "to groups via argument `K`.")
    }
    if (method == "ilp") {
      stop("The argument `K` should indicate the number of groups when using `method = 'ilp'`.")
    }
    if (argument_exists(categories)) {
      if (length(categories) != length(K)) {
        stop("Length of arguments `categories` and `K` differ, but ",
             "should be of same length (K can also be of length 1)")
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

  validate_input(method, "method", len = 1,
                 input_set = c("ilp", "exchange", "heuristic", "centroid"), not_na = TRUE)

  if (method == "ilp") {
    check_if_solver_is_available()
  }

  if (!inherits(objective, "function")) {
    validate_input(objective, "objective", input_set = c("distance", "diversity", "variance"), 
                   len = 1, not_na = TRUE)
    if (objective == "variance" && method == "ilp") {
      stop("You cannot use integer linear programming method to maximize the variance criterion. ",
           "Use objective = 'distance', or method = 'exchange' instead")
    }
    if (objective == "variance" && is_distance_matrix(x)) {
        stop("You cannot use a distance matrix with the objective 'variance'.")
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
#'
#' @return NULL
#'
#' @noRd

validate_input <- function(obj, argument_name, len = NULL, greater_than = NULL,
                           must_be_integer = FALSE, groupsize = NULL, input_set = NULL, 
                           objmode = NULL, not_na = FALSE) {

  self_validation(argument_name, len, greater_than,
                  must_be_integer, groupsize, input_set,
                  objmode, not_na)

  argument_name <- paste("Argument", argument_name)

  ## - Check length of input
  if (argument_exists(len)) {
    if (NROW(obj) != len) {
      stop(argument_name, " must have length ", len)
    }
  }

  ## - Check if input has to be greater than some value
  if (argument_exists(greater_than)) {
    if (any(obj <= greater_than)) {
      stop(argument_name, " must be greater than ", greater_than)
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

  ## - Check mode of input
  if (argument_exists(objmode)) {
    if (mode(obj) != objmode) {
      stop(argument_name, " must be ", objmode, ", but is ", mode(obj))
    }
  }

  return(invisible(NULL))
}

## Validate feature input
validate_data_matrix <- function(x) {
  x <- as.matrix(x)
  if (mode(x) != "numeric") {
    stop("Your data (the first argument `x`) should only contain numeric entries, but this is not the case.")
  }
  if (any(is.na(x))) {
    stop("Your data contains `NA`. I cannot proceed because ",
         "I cannot estimate similarity for data that has missing values. Sorry!")
  }
}

## Validate input for the `validate_input` function (these errors are
## not for users, but only for developers)
self_validation <- function(argument_name, len, greater_than,
                            must_be_integer, groupsize,
                            input_set, objmode, not_na) {

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
  validate_data_matrix(x)
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
  if (!requireNamespace("Rglpk", quietly = TRUE)) {
    stop("\n\nAn exact method was requested, but the package 'Rglpk' is not \n",
         "available. Possibly, you have to install the GNU linear \nprogramming kit first: \n\n",
         "- On windows, visit ",
         "http://gnuwin32.sourceforge.net/packages/glpk.htm \n\n",
         "- Use homebrew to install it on mac, 'brew install glpk' \n\n",
         "- 'sudo apt install libglpk-dev' on Ubuntu ",
         "\n\nThen, install the Rglpk package via ",
         "`install.packages('Rglpk')`.")
  }
  return(invisible(NULL))
}
