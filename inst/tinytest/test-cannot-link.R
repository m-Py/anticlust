
# Test if cannot_link constraint works
a <- anticlustering(rnorm(6), K = 2, cannot_link = cbind(1, 2))
expect_true(a[1] != a[2])

a <- anticlustering(rnorm(6), K = 2, cannot_link = cbind(1, 2), method = "ilp")
expect_true(a[1] != a[2])

a <- anticlustering(rnorm(6), K = 2, cannot_link = cbind(1, 2), method = "brusco")
expect_true(a[1] != a[2])


# now make it harder: find instance where two items are linked in optimal solution, 
# and then insert cannot-link constraint!

N <- 14
M <- 3
data <- matrix(rnorm(N*M), ncol = M)
a <- anticlustering(data, K = 2, method = "ilp")
one <- which(a == 1)[1]
other <- which(a == 1)[2]
diversity_objective(data, a)

b <- anticlustering(data, K = 2, cannot_link = cbind(one, other))
expect_true(b[one] != b[other])
expect_true(diversity_objective(data, a) >= diversity_objective(data, b))

b <- anticlustering(data, K = 2, cannot_link = cbind(one, other), method = "ilp")
expect_true(b[one] != b[other])
expect_true(diversity_objective(data, a) >= diversity_objective(data, b))

b <- anticlustering(data, K = 2, cannot_link = cbind(one, other), method = "brusco")
expect_true(b[one] != b[other])
expect_true(a[one] == a[other]) # this is necessary, but maybe we appreciate the reminder.
expect_true(diversity_objective(data, a) >= diversity_objective(data, b))



# Now: use multiple repetitions!

N <- 14
M <- 3
data <- matrix(rnorm(N*M), ncol = M)
a <- anticlustering(data, K = 2, method = "ilp")
one <- which(a == 1)[1]
other <- which(a == 1)[2]
diversity_objective(data, a)

b <- anticlustering(data, K = 2, cannot_link = cbind(one, other), repetitions = 10)
expect_true(b[one] != b[other])
expect_true(diversity_objective(data, a) >= diversity_objective(data, b))

b <- anticlustering(data, K = 2, cannot_link = cbind(one, other), method = "brusco", repetitions = 10)
expect_true(b[one] != b[other])
expect_true(a[one] == a[other]) # this is necessary, but maybe we appreciate the reminder.
expect_true(diversity_objective(data, a) >= diversity_objective(data, b))



## TEST FAILS: Use method to fill elements from must-link branch!
N <- 140
M <- 5
data <- matrix(rnorm(N*M), ncol = M)
b <- anticlustering(data, K = 2, cannot_link = cbind(1:2, 2:3), repetitions = 10, method = "local-maximum")
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])
b <- anticlustering(data, K = 2, cannot_link = cbind(1:2, 2:3), method = "brusco", repetitions = 10)
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])

