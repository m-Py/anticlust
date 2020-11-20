

N <- 15
M <- 2
K <- 3
data <- matrix(rnorm(N * M), ncol = M)
tt <- Sys.time()
foo_clust(data, K)
Sys.time() - tt
dists <- distances_from_centroid(data)
tt <- Sys.time()
nns <- nn2(data, k = N / K)
Sys.time() - tt
nns$nn.idx[which.min(dists), ]
