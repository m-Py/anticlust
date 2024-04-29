
#test setup of ILP
library("anticlust")

context("ILP setup")

#anticlustering_ilp
test_that("ILP is set up as expected", {
  conditions <- expand.grid(m = 1:4, p = 2:4)
  for (k in 1:nrow(conditions)) {
    m_features <- conditions[k, "m"]
    p_anticlusters <- conditions[k, "p"]
    n_elements <- p_anticlusters * 2 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    distances <- as.matrix(dist(features))
    ilp <- anticlustering_ilp(distances, p_anticlusters)
    expect_equal(nrow(ilp$constraints), choose(n_elements, 3) * 3 + n_elements)
    expect_equal(sum(is.na(ilp$constraints)), 0)

    ## Test that distance matrix is correctly transfered into objective vector
    distances <- as.matrix(distances)
    costs <- ilp$costs
    for (i in 1:n_elements) {
      for (j in 1:n_elements) {
        if (i >= j) next
        expect_equal(costs[costs$i == i & costs$j == j, ]$costs, distances[i, j])
      }
    }
  }
})
