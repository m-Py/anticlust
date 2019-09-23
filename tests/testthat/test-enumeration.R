
context("Test complete enumeration implementation")
library("anticlust")

test_that("complete enumeration computes optimal solution", {
  conditions <- expand.grid(m = 1:4, p = 2:3)
  for (i in nrow(conditions)) {
    m_features <- conditions[i, "m"]
    K <- conditions[i, "p"]
    n_elements <- K * 4 # n must be multiplier of p
    features <- matrix(rnorm(n_elements * m_features), ncol = m_features)
    enum_anticlusters <- enum_anticlustering(features, K = K)
    # Use preclustering as resticting information in anticlustering
    ilp_anticlusters <- anticlustering(features, K = K, preclustering = FALSE,
                                       method = "ilp")
    enum_obj <- round(obj_value_distance(enum_anticlusters, features), 10)
    ilp_obj  <- round(obj_value_distance(ilp_anticlusters, features), 10)
    expect_equal(enum_obj, ilp_obj)
  }
})
