
library("anticlust")
# Preclustering works with all criteria

features <- schaper2019[, 3:6]
matches <- matching(features, p = 3)

for (obj in c("variance", "kplus", "dispersion", "diversity")) {
  anticlusters <- anticlustering(
    features, 
    K = 3, 
    objective = obj,
    preclustering = TRUE
  )
  
  expect_true(all(table(matches, anticlusters) == 1))
}

# kplus_anticlustering must work as well 
anticlusters <- kplus_anticlustering(features, K = 3, preclustering = TRUE)

expect_true(all(table(matches, anticlusters) == 1))
