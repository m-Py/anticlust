

N <- 20
M <- 2
data <- matrix(rnorm(N * M), ncol = M)
foo_clust(data, 1)
colMeans(data)
data[order(distances_from_centroid(data)), ]
