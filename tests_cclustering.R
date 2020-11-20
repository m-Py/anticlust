

N <- 10
M <- 2
K <- 2
data <- matrix(rnorm(N * M), ncol = M)
foo_clust(data, K)
dists <- distances_from_centroid(data)
nn2(data, k = 5)$nn.idx[which.min(dists), ]
