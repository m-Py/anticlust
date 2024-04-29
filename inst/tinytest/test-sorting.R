

library("anticlust")

# data table sorting works
mat <- cbind(10:1, 1:10)
colnames(mat) <- c("one", "two")
expect_equal(mat[, 1], 10:1)
expect_equal(mat[, 2], 1:10)
expect_equal(anticlust:::sort_by_col(mat, 1)[, 1], 1:10)
expect_equal(anticlust:::sort_by_col(mat, 1)[, "one"], 1:10)
expect_equal(anticlust:::sort_by_col(mat, 1)[, 2], 10:1)
expect_equal(anticlust:::sort_by_col(mat, 1)[, "two"], 10:1)
expect_error(anticlust:::sort_by_col(mat, 1:2))
expect_error(anticlust:::sort_by_col(mat, c(TRUE, FALSE)))
expect_error(anticlust:::sort_by_col(mat, c("TRUE", "FALSE")))
expect_error(anticlust:::sort_by_col(1, 1))
