
#' Get binary representation of categorical variables
#'
#' @param categories A vector, data.frame or matrix representing one
#'     or several categorical variables
#' @param use_combinations Logical, should the output also include columns representing
#'    the combination / interaction of the categories (defaults to \code{FALSE}).
#'
#' @return A matrix encoding the categorical variable(s) in binary form. 
#'
#' @details
#' 
#' The conversion of categorical variables to binary variables is done via
#' \code{\link[stats]{model.matrix}}. Since version 0.8.9, each category
#' of a categorical variable is coded by a separate variable. So this is not
#' 'dummy' coding, which is often used to encode predictors in statistical 
#' analysis. Dummy coding uses a reference category that has only zeros for 
#' each variable, while all other categories consist of a 1 and otherwise zeros. 
#' This implies that there is a different distance to the reference category 
#' than among the other categories, which is unwarranted in anticlustering.
#' 
#' This function can be used to include categorical variables as part of the 
#' optimization criterion in anticlustering, rather than including them as hard constraints as done when using the 
#' argument \code{categories} in \code{\link{anticlustering}} (or \code{\link{fast_anticlustering}}). 
#' This way, categorical variables are treated as numeric variables, 
#' which can be useful when there are several
#' categorical variables or when the group sizes are unequal (or both).
#' See examples. Please see the vignette 'Using categorical variables with anticlustering'
#' for more information on this approach.
#' 
#' @importFrom stats as.formula model.matrix contrasts
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
#' @export
#'
#' @examples
#' 
#' # How to encode a categorical variable with three levels:
#' unique(iris$Species)
#' categories_to_binary(iris$Species)[c(1, 51, 101), ]
#' 
#' # Use Schaper data set for anticlustering example
#' data(schaper2019)
#' features <- schaper2019[, 3:6]
#' K <- 3
#' N <- nrow(features) 
#' 
#' # - Generate data input for k-means anticlustering -
#' # We conduct k-plus anticlustering by first generating k-plus variables, 
#' # and also include the categorical variable as "numeric" input for the 
#' # k-means optimization (rather than as input for the argument \code{categories})
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
  model.matrix(
    as.formula(formula_string), 
    data = categories,
    contrasts.arg = lapply(categories, contrasts, contrasts=FALSE) # this ensures that each level of the category has a binary variable
  )[ ,-1, drop = FALSE]
}
