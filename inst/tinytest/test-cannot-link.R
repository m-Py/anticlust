
library(anticlust)
library(tinytest)

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


## Test for larger data sets:
N <- 140
M <- 5
data <- matrix(rnorm(N*M), ncol = M)
b <- anticlustering(data, K = 2, cannot_link = cbind(1:2, 2:3), repetitions = 10, method = "local-maximum")
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])
b <- anticlustering(data, K = 2, cannot_link = rbind(1:2, 2:3), method = "brusco", repetitions = 10, objective = "variance")
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])

expect_true(table(b)[1] == table(b)[2])

# different sized groups:
b <- anticlustering(data, K = c(100, 40), cannot_link = rbind(1:2, 2:3), method = "brusco", repetitions = 10)
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])

expect_true(sort(table(b))[1] == 40)
expect_true(sort(table(b))[2] == 100)

# more constraints, different sized groups
b <- anticlustering(data, K = c(100, 40), cannot_link = rbind(1:2, 2:3, 4:5, 6:7), method = "brusco", repetitions = 10)
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])
expect_true(b[4] != b[5])
expect_true(b[6] != b[7])

expect_true(sort(table(b))[1] == 40)
expect_true(sort(table(b))[2] == 100)

# more groups, different sized groups
groups <- c(20, 20, 40, 30, 10, 10, 10)
b <- anticlustering(data, K = groups, cannot_link = rbind(1:2, 2:3, 4:5, 6:7), method = "brusco", repetitions = 10)
expect_true(b[1] != b[2])
expect_true(b[2] != b[3])
expect_true(b[4] != b[5])
expect_true(b[6] != b[7])

expect_true(all(sort(groups) == sort(table(b))))

# Test that all argument combinations go through without error

N <- 30
M <- 5
K <- 5
data <- matrix(rnorm(N*M), ncol = M)

combinations <- expand.grid(
  objective = c("kplus", "variance", "diversity", "average-diversity"),
  repetitions = c(1, 10),
  method = c("brusco", "local-maximum", "exchange"),
  standardize = c(TRUE, FALSE)
)


for (i in 1:nrow(combinations)) {
  objective <- combinations[i, "objective"]
  method <- combinations[i, "method"]
  repetitions <- combinations[i, "repetitions"]
  standardize <- combinations[i, "standardize"]
  a <- anticlustering(data, K = K, objective = objective, method = method, repetitions = repetitions, standardize = standardize) # Once without cannot-link constraints
  b <- anticlustering(data, K = K, objective = objective, method = method, repetitions = repetitions, standardize = standardize, cannot_link = rbind(1:2, 2:3))
  expect_true(b[1] != b[2])
  expect_true(b[2] != b[3])
}

## Expect errors

expect_error(
  anticlustering(data, K = K, cannot_link = rbind(1:2, 2:3), categories = sample(LETTERS, size = N, replace = TRUE)),
  pattern = "categories"
)
expect_error(
  anticlustering(data, K = K, cannot_link = rbind(1:2, 2:3),preclustering = TRUE),
  pattern = "preclustering"
)

expect_error(
  anticlustering(data, K = K, cannot_link = rbind(1:3, 5:7)),
  pattern = "columns"
)

# Test cannot-link constraints via vector

cannot_link <- sample(N/2, size = N, replace = TRUE)
K <- max(table(cannot_link))
a <- anticlustering(data, K = K, cannot_link = cannot_link)

expect_true(all(table(a, cannot_link) <= 1))

expect_error(anticlustering(data, K = K-1, cannot_link = cannot_link))
