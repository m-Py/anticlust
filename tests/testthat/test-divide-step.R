

# Test 1: per vector
x <- sample(5:9, size = 20, replace = TRUE)
u <- categorize_vector(x, 3)

listed <- tapply(x, u, c)
expect_true(max(listed[[1]]) <= min(listed[[2]]))
expect_true(max(listed[[2]]) <= min(listed[[3]]))


# Test 2: for a data frame
mat <- matrix(sample(5:9, size = 20, replace = TRUE), ncol = 2)
splitted <- split_data(mat, c(2, 2))

sort_by_col(sort_by_col(cbind(mat, splitted), 3), 4)
