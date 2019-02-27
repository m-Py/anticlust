
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
#'   refer to the "mode" of the matrix. A special case is the
#'   \code{class_string} "groupvar" that is expanded to
#'   c("factor", "character", "numeric", "integer", "double").
#' @param len Optional numeric vector for objects having a length
#'   (mostly for vectors).
#' @param gt0 Optional logical vector indicating if numeric input has
#'   to be greater than 0.
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
                           len = NULL, gt0 = FALSE, must_be_integer = FALSE, groupsize = NULL,
                           input_set = NULL, objmode = NULL, not_na = FALSE) {

  self_validation(argument_name, class_string, len, gt0,
                  must_be_integer, groupsize, input_set, objmode)

  ## - Check class of object
  if (argument_exists(class_string))  {
    # Allow for all numeric types:
    if ("numeric" %in% class_string) {
      class_string <- c(class_string, "integer", "double")
    }
    # Case - grouping variable: Allow for numeric, character or factor
    if ("groupvariable" %in% class_string) {
      class_string <- setdiff(c(class_string, "factor", "character",
                                "numeric", "integer", "double"),
                              "groupvariable")
    }
    correct_class <- class(obj) %in% class_string
    if (!correct_class) {
      stop(argument_name, " must be of class '",
           paste(class_string, collapse = "' or '"), "'")
    }
  }

  ## - Check length of input
  if (argument_exists(len)) {
    if (length(obj) != len) {
      stop(argument_name, " must have length ", len)
    }
  }

  ## - Check if input has to be greater than 0
  if (gt0 == TRUE && any(obj <= 0)) {
    stop(argument_name, " must be greater than 0")
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
      stop(argument_name, " must be integer of mode ", objmode,
           ", but is of mode ", mode(obj))
    }
  }

  if (not_na == TRUE) {
    if (is.na(obj)) {
      stop(argument_name, " cannot must not be NA but is NA")
    }
  }

  return(invisible(NULL))
}

## Validate input for the `validate_input` function (these errors are
## not for users, but only for developers)
self_validation <- function(argument_name, class_string, len, gt0,
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
  stopifnot(gt0 %in% c(TRUE, FALSE))
  stopifnot(length(gt0) == 1)
  stopifnot(must_be_integer %in% c(TRUE, FALSE))
  stopifnot(length(must_be_integer) == 1)
  stopifnot(not_na %in% c(TRUE, FALSE))
  stopifnot(length(not_na) == 1)

  if (argument_exists(groupsize)) {
    stopifnot(mode(groupsize) == "numeric")
    stopifnot(length(groupsize) == 1)
  }

  if (argument_exists(input_set) && len != 1) {
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
# anticlusters. I think this function is only used in tests.
# maybe replace / delete it later.
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
