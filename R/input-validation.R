

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
  N <- nrow(as.matrix(x))

  validate_input(preclustering, "preclustering", "logical", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE)

  if (argument_exists(categories) && preclustering == TRUE) {
    stop("It is currently not possible to apply both preclustering and categorical constraints.")
  }

  # Allow that K is an initial assignment of elements to clusters
  if (length(K) == 1) {
    validate_input(K, "K", len = 1, greater_than = 1, must_be_integer = TRUE, not_na = TRUE)
  } else {
    validate_input(K, "K", "numeric", len = N, not_na = TRUE)
    if (method != "exchange") {
      stop("Passing an initial cluster assignment via the argument `K` ",
           "only works with method = 'exchange'")
    }
    if (argument_exists(categories)) {
      if (length(categories) != length(K)) {
        stop("Length of arguments `categories` and `K` differ, but should be of same length (K can also be of length 1)")
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
    if (preclustering == TRUE) {
      stop("The data set cannot be clustered into equal-sized parts.")
    }
  }

  validate_input(method, "method", len = 1,
                 input_set = c("ilp", "exchange", "heuristic"), not_na = TRUE)

  if (method == "ilp") {
    solver <- solver_available()
    if (solver == FALSE) {
      stop("\n\nAn exact solution was requested, but none of the linear ",
           "programming \npackages 'Rglpk', 'gurobi', or 'Rcplex' is ",
           "available. \n\nTry `method = 'exchange'` or install ",
           "a linear programming solver \nto obtain an exact solution. ",
           "For example, install the GNU linear \nprogramming kit: \n\n",
           "- On windows, visit ",
           "http://gnuwin32.sourceforge.net/packages/glpk.htm \n\n",
           "- Use homebrew to install it on mac, 'brew install glpk' \n\n",
           "- 'sudo apt install libglpk-dev' on Ubuntu ",
           "\n\nThen, install the Rglpk package via ",
           "`install.packages('Rglpk')`. \n\nOtherwise, you may obtain ",
           "a license for one of ",
           "the commercial solvers \ngurobi or IBM CPLEX (they are free ",
           "for academic use).")
    }
  }

  if (class(objective) != "function") {
    validate_input(objective, "objective", input_set = c("distance", "variance"), len = 1, not_na = TRUE)
    if (objective == "variance" && method == "ilp") {
      stop("You cannot use integer linear programming method to maximize the variance criterion. ",
           "Use objective = 'distance', or method = 'exchange' instead")
      if (objective == "variance" && is_distance_matrix(x)) {
        stop("You cannot use a distance matrix with the objective 'variance'.")
      }
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
#' @param class_string A character vector of legal classes. If
#'   \code{class_string} is "numeric", it will be expanded to
#'   c("numeric", "integer", "double"). The class is tested via the
#'   function \code{class}. This means that if \code{obj} is a matrix,
#'   it is necessary to pass \code{class_string = "matrix"}; you cannot
#'   refer to the "mode" of the matrix.
#' @param len Optional numeric vector for objects having a length
#'   (mostly for vectors).
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

validate_input <- function(obj, argument_name, class_string = NULL,
                           len = NULL, greater_than = NULL, must_be_integer = FALSE,
                           groupsize = NULL, input_set = NULL, objmode = NULL,
                           not_na = FALSE) {

  self_validation(argument_name, class_string, len, greater_than,
                  must_be_integer, groupsize, input_set,
                  objmode, not_na)

  argument_name <- paste("Argument", argument_name)
  
  ## - Check class of object
  if (argument_exists(class_string))  {
    # Allow for all numeric types:
    if ("numeric" %in% class_string) {
      class_string <- c(class_string, "integer", "double")
    }
    correct_class <- class(obj) %in% class_string
    if (!correct_class) {
      stop(argument_name, " must be of class '",
           paste(class_string, collapse = "' or '"), "'")
    }
  }

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
      stop(argument_name, " must be ", objmode,
           ", but is ", mode(obj))
    }
  }

  return(invisible(NULL))
}

## Validate feature input
validate_data_matrix <- function(x) {
  x <- as.matrix(x)
  validate_input(x, "features", objmode = "numeric")
  if (sum(!complete.cases(x)) >= 1) {
    warning("There are missing values in your data, take care!")
  }
}

## Validate input for the `validate_input` function (these errors are
## not for users, but only for developers)
self_validation <- function(argument_name, class_string, len, greater_than,
                            must_be_integer, groupsize,
                            input_set, objmode, not_na) {
  if (argument_exists(class_string)) {
    stopifnot(class(class_string) == "character")
    stopifnot(class(argument_name) == "character")
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
