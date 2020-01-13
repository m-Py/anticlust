
### Input validation for stimulus selection functions
validate_input_selection <- function(data, split_by, equalize, balance, design, n) {
  if (argument_exists(split_by)) {
    validate_input(design, "design", must_be_integer = TRUE, not_na = TRUE)
    if (length(design) != length(split_by)) {
      stop("Length of argument `design` must match length of argument `split_by`.")
    }
  } else {
    validate_input(design, "design", len = 1, must_be_integer = TRUE, not_na = TRUE)
  }
  
  if (argument_exists(n)) {
    validate_input(n, "n", len = 1, must_be_integer = TRUE, not_na = TRUE)
  }
  
  # test that legal colnames are given
  cols <- c(split_by, equalize, balance)
  good_colnames <- all(cols %in% colnames(data))
  if (good_colnames == FALSE) {
    stop("The arguments `split_by`, `equalize` and `balance` must refer to column names in your data.")
  }
  return(invisible(NULL))
}
