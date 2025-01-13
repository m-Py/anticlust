
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

## Third clique has size of 3. 
# THIS IS CURRENTLY A FAILING TEST. ENSURE THAT THE ALGORITHM DOES NOT ALWAYS GENERATE THE SAME
# COMBINATION OF EXCHANGE PARTNERS. IT MATCHES THE FIRST FIT (which is the first clique of size 2 + the first singleton on position 8)
# INSTEAD OF A RANDOM FIT. 
anticlust:::get_exchange_partners_clique(cliques, 3, init, must_link)[["sample_ids"]]
expect_true(
  any(!replicate(1000, anticlust:::get_exchange_partners_clique(cliques, 3, init, must_link)[["sample_ids"]]) == c(1, 2, 8))
)

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

tt <- anticlustering(1:10, K = 2, must_link = must_link)
tt2 <- anticlustering(1:10, K = 2, must_link = must_link, method = "2PML")

expect_true(must_link_constraints_valid(tt, must_link))
expect_true(must_link_constraints_valid(tt2, must_link))

# test constructing merged clusters from must-link constraints and its reversal
expect_true(all(
  tt == anticlust:::merged_cluster_to_original_cluster(
    anticlust:::original_cluster_to_merged_cluster(tt, must_link), 
    must_link
  )
))


