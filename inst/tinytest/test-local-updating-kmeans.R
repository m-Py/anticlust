
library("anticlust")

# local updating works correctly for k-means anticlustering

M <- 2
N <- 60
K <- 2
clusters <- sample(rep(1:K, c(40, 20)))
features <- matrix(rnorm(N * M), ncol = M)

# Swap two items and check objective
to_swap <- sample(1:K, size = 2)
swap1 <- which(clusters == to_swap[1])[1]
swap2 <- which(clusters == to_swap[2])[1]

centers <- anticlust:::cluster_centers(features, clusters)

# Compute clusters after swap
local_update_centers <- anticlust:::update_centers(
  centers, 
  features, 
  swap1, 
  swap2, 
  clusters[swap1], 
  clusters[swap2], 
  table(clusters)
)

# swap clusters, and recompute cluster centers
tmp <- clusters[swap1]
clusters[swap1] <- clusters[swap2]
clusters[swap2] <- tmp

centers_after_swap <- anticlust:::cluster_centers(features, clusters)

expect_equal(local_update_centers, centers_after_swap)
