

N <- 1000
M <- 2
K <- 20
data <- matrix(rnorm(N * M), ncol = M)
tt <- Sys.time()
this <- foo_clust(data, K)
Sys.time() - tt
tt <- Sys.time()
nns <- balanced_clustering(data, K)
Sys.time() - tt
#nns$nn.idx[which.min(dists), ]
par(mfrow = c(1, 2))
plot_clusters(data, clusters = this)
plot_clusters(data, clusters = nns)
diversity_objective(data, this)
diversity_objective(data, nns)

# Cluster a data set and visualize results
N <- 1000
lds <- data.frame(f1 = rnorm(N), f2 = rnorm(N))
cl <- foo_clust(lds, K = 10)
cl2 <- balanced_clustering(lds, K = 10)
par(mfrow = c(1, 2))
plot_clusters(lds, clusters = cl)
plot_clusters(lds, clusters = cl2)
