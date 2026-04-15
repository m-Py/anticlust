
library(anticlust)
library(tinytest)

# test that one-hot encoding is done correctly

generate_categorical_data <- function(N, M, P, prob = NULL) {
  data <- data.frame(matrix(sample(P, replace = TRUE, size = N*M, prob = prob), ncol = M))
  # wandel Datensatz in Faktor um (jede Spalte)
  as.data.frame(lapply(data, as.factor))
}

# old version of categories_to_binary(), that clearly uses one-hot encoding:
one_hot_encoding <- function(categories, use_combinations = FALSE) {
  categories <- data.frame(categories)
  combine_by <- ifelse(use_combinations, " * ", " + ")
  formula_string <- paste("~", paste(colnames(categories), collapse = combine_by), collapse = "")
  model.matrix(
    as.formula(formula_string),
    data = categories,
    contrasts.arg = lapply(categories, contrasts, contrasts=FALSE) # this ensures that each level of the category has a binary variable
  )[ ,-1, drop = FALSE]
}


d <- generate_categorical_data(20, 4, 5)

one_hot1 <- one_hot_encoding(d)
dim(one_hot1)
one_hot2 <- anticlust::categories_to_binary(d)
dim(one_hot2)

expect_equal(one_hot1, one_hot2)
