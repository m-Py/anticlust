
#' Get binary representation of categorical variables
#'
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical variables
#' @param use_combinations Logical, should the output also include columns representing
#'    the combination / interaction of the categories (defaults to \code{FALSE}).
#'
#' @return A matrix representing the categorical variables in binary form ("dummy coding")
#'
#' @details
#' 
#' The conversion of categorical variable to binary variables is done via
#' \code{\link[stats]{model.matrix}}. This function can be used to include 
#' categorical variables as part of the optimization criterion in k-means / 
#' k-plus anticlustering, rather than including them as hard constraints as 
#' done in \code{\link{anticlustering}}. This can be useful when there are several
#' categorical variables or when the group sizes are unequal (or both).
#' See examples.
#' 
#' @importFrom stats as.formula
#'
#'
#' @author
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#'
#' @references
#' 
#' Papenberg, M. (2024). K-plus Anticlustering: An Improved k-means Criterion for 
#' Maximizing Between-Group Similarity. British Journal of Mathematical and 
#' Statistical Psychology, 77(1), 80--102. https://doi.org/10.1111/bmsp.12315
#' 
#' 
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' 
#' # Use Schaper data set for example
#' data(schaper2019)
#' features <- schaper2019[, 3:6]
#' K <- 3
#' N <- nrow(features) 
#' 
#' # - Generate data input for k-means anticlustering -
#' # We conduct k-plus anticlustering by first generating k-plus variables, 
#' # and also include the categorical variable as "numeric" input for the 
#' # k-means optimization (rather than as input for the argument `categories`)
#' 
#' input_data <- cbind(
#'   kplus_moment_variables(features, T = 2), 
#'   categories_to_binary(schaper2019$room) 
#' )
#' 
#' kplus_groups <- anticlustering(
#'   input_data, 
#'   K = K,
#'   objective = "variance",
#'   method = "local-maximum", 
#'   repetitions = 10
#' )
#' mean_sd_tab(features, kplus_groups)
#' table(kplus_groups, schaper2019$room) # argument categories was not used!
#' 
#'  

categories_to_binary <- function(categories, use_combinations = FALSE) {
  validate_input(use_combinations, "use_combinations", objmode = "logical", len = 1,
                 input_set = c(TRUE, FALSE), not_na = TRUE, not_function = TRUE)
  categories <- data.frame(categories)
  categories <- as.data.frame(lapply(categories, as.factor))
  combine_by <- ifelse(use_combinations, " * ", " + ")
  formula_string <- paste("~", paste(colnames(categories), collapse = combine_by), collapse = "")
  model.matrix(as.formula(formula_string), data = categories)[ ,-1, drop = FALSE]
}
