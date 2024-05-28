
# Test if cannot_link constraint works
a <- anticlustering(rnorm(6), K = 2, cannot_link = cbind(1, 2))
expect_true(a[1] != a[2])

a <- anticlustering(rnorm(6), K = 2, cannot_link = cbind(1, 2), method = "ilp")
expect_true(a[1] != a[2])

a <- anticlustering(rnorm(6), K = 2, cannot_link = cbind(1, 2), method = "brusco")
expect_true(a[1] != a[2])
