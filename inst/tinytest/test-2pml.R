
library("anticlust")
library("tinytest")

must_link <- c(1, 1, 2, 2, 3, 3, 3, 4, 5, 6)

ml_indices <- anticlust:::get_must_link_indices(must_link)
singletons <- anticlust:::get_singletons(ml_indices)
cliques <- anticlust:::get_cliques(ml_indices)
expect_true(all(singletons == 8:10))
expect_true(all(cliques[[1]] == 1:2))
expect_true(all(cliques[[2]] == 3:4))
expect_true(all(cliques[[3]] == 5:7))

## Test generating exchange partners

init <- anticlust:::merged_cluster_to_original_cluster(
  anticlust:::init_must_link_groups(
    length(must_link), 
    IDs_initial = must_link, 
    IDs_reduced = anticlust:::get_must_link_indices(must_link), 
    target_groups = c(5, 5)
  ), must_link
)

# the last 3 elements can serve as exchange partners for a clique of size 3
ep <- anticlust:::get_exchange_partners_singletons(singletons, init[singletons], 3)
expect_true(all(sort(ep$sample_ids) == 8:10))

# now test for clique of size 2 (any combination of the last 3 works)
ep <- anticlust:::get_exchange_partners_singletons(singletons, init[singletons], 2)
expect_true(all(ep$sample_ids %in% 8:10))

## test generating exchange partners from cliques (and singletons)
# (this does not work because clique must be larger than 2)
expect_error(anticlust:::get_exchange_partners_clique(cliques, 1, init, must_link))

## Third clique has size of 3. It must be matched against positions 1 and 2, and (8, 9, or 10)

anticlust:::get_exchange_partners_clique(cliques, 3, init, must_link)[["sample_ids"]]
tt <- replicate(1000, anticlust:::get_exchange_partners_clique(cliques, 3, init, must_link)[["sample_ids"]])

# We have two cliques of size 2. (1) Indices 1, 2; (2) Indices 3, 4.
# These two cliques must be partnered with (a) the clique of size 3, (b) the 3 singletons.

# Test that only those matches occur are matched
expect_true(all(apply(tt, 2, function(x) x == c(1:2, 8) | x == c(1:2, 9) | x == c(1:2, 10) | x == c(3:4, 8) | x == c(3:4, 9) | x == c(3:4, 10))))

## But not always c(1:2, 8) or c(3:4, 8)! (which we get when using `base::match()` rather than `anticlust:::random_match()`)
expect_true(any(!tt == c(c(1:2, 8))))
expect_true(any(!tt == c(c(3:4, 8))))

# Test random match function. The match should always include 5 and 1, but 2, 3 or 4 is arbitrary.
combination <- c(1, 2, 3)
frequencies <- c(2, 1, 1, 1, 3)
match(combination, frequencies) # always c(2, 1, 5); we also like to have c(3, 1, 5), c(4, 1, 5)
# However, positions 1 and 5 must always be matched
expect_true(all(replicate(100, 1 %in% anticlust:::random_match(combination, frequencies))))
expect_false(all(replicate(100, 2 %in% anticlust:::random_match(combination, frequencies))))
expect_true(all(replicate(100, 5 %in% anticlust:::random_match(combination, frequencies))))


## Test that a 4 clique can be exchanged with 2, 2 ### FAILS!

must_link <- c(1, 1, 1, 1, 5, 2, 2, 3, 3, 4)

ml_indices <- anticlust:::get_must_link_indices(must_link)
cliques <- anticlust:::get_cliques(ml_indices)
init <- anticlust:::merged_cluster_to_original_cluster(
  anticlust:::init_must_link_groups(
    length(must_link), 
    IDs_initial = must_link, 
    IDs_reduced = anticlust:::get_must_link_indices(must_link), 
    target_groups = c(5, 5)
  ), must_link
)

exchanges <- anticlust:::get_exchange_partners_clique(cliques, 1, init, must_link)
expect_true(!is.null(exchanges$cluster_id))
expect_true(all(exchanges$sample_ids == 6:9))

## Test that a 5 clique can be exchanged with 2, 2, 1

must_link <- c(1, 1, 1, 1, 1, 2, 2, 3, 3, 4)

ml_indices <- anticlust:::get_must_link_indices(must_link)
cliques <- anticlust:::get_cliques(ml_indices)
init <- anticlust:::merged_cluster_to_original_cluster(
  anticlust:::init_must_link_groups(
    length(must_link), 
    IDs_initial = must_link, 
    IDs_reduced = anticlust:::get_must_link_indices(must_link), 
    target_groups = c(5, 5)
  ), must_link
)

exchanges <- anticlust:::get_exchange_partners_clique(cliques, 1, init, must_link)
expect_true(!is.null(exchanges$cluster_id))
expect_true(all(exchanges$sample_ids == 6:10))

## Test stuff

N <- 100
must_link <- sample(N, N, replace = TRUE)

ml_indices <- anticlust:::get_must_link_indices(must_link)
cliques <- anticlust:::get_cliques(ml_indices)
init <- anticlust:::merged_cluster_to_original_cluster(
  anticlust:::init_must_link_groups(
    length(must_link), 
    IDs_initial = must_link, 
    IDs_reduced = anticlust:::get_must_link_indices(must_link), 
    target_groups = c(N/2, N/2)
  ), must_link
)

exchanges <- anticlust:::get_exchange_partners_clique(cliques, which(sapply(cliques, function(x) length(x) > 2))[1], init, must_link)
# there must at least be one clique among the exchange partners:
expect_true(any(duplicated(must_link[exchanges$sample_ids])))


########################################
### Now do some general testing on the must-link constraints. E.g., are they valid after anticlustering?

# Function that test if constraints are valid
must_link_constraints_valid <- function(cl, must_link) {
  same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
  all_good <- c()
  for (i in same) {
    all_good <- c(all_good, all(cl[must_link == i] == cl[must_link == i][1]))
  }
  all(all_good)
}

N <- 10

must_link <- sample(N, size = N, replace = TRUE)

tt <- anticlustering(1:N, K = 2, must_link = must_link)
tt2 <- anticlustering(1:N, K = 2, must_link = must_link, method = "2PML")

expect_true(must_link_constraints_valid(tt, must_link))
expect_true(must_link_constraints_valid(tt2, must_link))

# test constructing merged clusters from must-link constraints and its reversal
expect_true(all(
  tt == anticlust:::merged_cluster_to_original_cluster(
    anticlust:::original_cluster_to_merged_cluster(tt, must_link), 
    must_link
  )
))



## now a random larger data set

K <- 20
group_sizes <- 12
N <- group_sizes * K # here 128
M <- 5
data <- matrix(rnorm(N * M), ncol = M)
distances <- as.matrix(dist(data))

# Generate random patient IDs
must_link <- sample(N, replace = TRUE)


tt1 <- tryCatch(
  anticlustering(distances, K, must_link = must_link, method = "2PML"),
  error = function(e) e
)

if (!"simpleError" %in% class(tt1)) {
  expect_true(must_link_constraints_valid(tt1, must_link))
}


tt2 <- tryCatch(
  anticlustering(distances, K, must_link = must_link, method = "2PML", repetitions = 40),
  error = function(e) e
)

if (!"simpleError" %in% class(tt2)) {
  expect_true(must_link_constraints_valid(tt2, must_link))
}

tt0 <- anticlustering(distances, K, must_link = must_link) # exchange method / Phase 1 of 2PML without ensured local optimality

if (!"simpleError" %in% class(tt0)) {
  diversity_objective(distances, tt0)
  diversity_objective(distances, tt1)
  diversity_objective(distances, tt2)
}
