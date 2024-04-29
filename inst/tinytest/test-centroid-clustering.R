

context("Centroid clustering")
library("anticlust")

# centroid clustering algorithm relies on rownames, this used to generate
# a warning for subsetted data. make sure this warning no longer occurs
test_that("centroid clustering does not throw warning", {
  # does it work when the elements have rownames?
  features <- rnorm(100)
  features <- matrix(features, ncol = 2)
  rownames(features) <- paste0("foo", 1:nrow(features))
  expect_warning(
    cl <- balanced_clustering(
      dist(features),
      K = 2
    ),
    regexp = NA
  )

  ## does it work with subsetted data?
  features <- schaper2019[schaper2019$room == "kitchen", 3:4]
  expect_warning(
    cl <- balanced_clustering(
      features,
      K = 2
    ),
    regexp = NA
  )

})

