

context("Centroid clustering")
library("anticlust")

# centroid clustering algorithm relies on rownames, this used to generate
# a warning for subsetted data. make sure this warning no longer occurs
test_that("centroid clustering does not throw warning", {
  expect_warning(
    balanced_clustering(
      distances = dist(schaper2019[schaper2019$room == "kitchen", 3:4]),
      K = 2
    ),
    regexp = NA
  )
})
