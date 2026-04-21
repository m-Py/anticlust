
library("anticlust")
library("tinytest")

cov_mat <- cov(iris[, -5])
lambda <- 0
x <- iris[1, -5]
y <- iris[2, -5]

md1 <- mahalanobis_distance(
  x, y,
  cov_mat = cov_mat,
  lambda = lambda
)

inv_cov_mat <- MASS::ginv(cov_mat+lambda*diag(nrow(cov_mat)))

maha_dist <- function(x, y, inv_cov_mat) {
  delta <- as.numeric(x - y)
  as.numeric(t(delta) %*% inv_cov_mat %*% delta)
}

md2 <- maha_dist(x, y, inv_cov_mat)

expect_equal(md1, md2)

md3 <- mahalanobis_distance(
  x, y,
  cov_mat = cov_mat,
  lambda = 2
)

inv_cov_mat <- MASS::ginv(cov_mat+2*diag(nrow(cov_mat)))
md4 <- maha_dist(x, y, inv_cov_mat)
expect_equal(md3, md4)
