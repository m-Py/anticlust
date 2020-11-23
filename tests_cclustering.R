

N <- 10000
M <- 2
K <- 5
data <- matrix(rnorm(N * M), ncol = M)
tt <- Sys.time()
foo_clust(data, N / K)
Sys.time() - tt
tt <- Sys.time()
nns <- balanced_clustering(data, N / K)
Sys.time() - tt
#nns$nn.idx[which.min(dists), ]


# Cluster a data set and visualize results
N <- 1000
lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
cl <- foo_clust(lds, K = 10)
cl2 <- balanced_clustering(lds, K = 10)
par(mfrow = c(1, 2))
plot_clusters(lds, clusters = cl)
plot_clusters(lds, clusters = cl2)
