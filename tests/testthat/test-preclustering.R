
library("anticlust")

context("Are preclusters balanced across anticlusters")

test_that("Preclustering works with all criteria", {
  
  features <- schaper2019[, 3:6]
  
  for (obj in c("variance", "kplus", "dispersion", "diversity")) {
    anticlusters <- anticlustering(
      features, 
      K = 3, 
      objective = obj,
      preclustering = TRUE
    )
    
    matches <- matching(features, p = 3)
    expect_true(all(table(matches, anticlusters) == 1))
  }
  
})
