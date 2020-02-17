# anticlust 0.4.1

## Minor 

- Some changes to documentation
- There is now a package website at https://m-py.github.io/anticlust/
- TODO: Additional error handling

## Internal

- Improved efficiency of k-means anticlustering: on each exchange 
iteration, only recomutes distances from clusters whose elements
have been swapped (mostly relevant for larger K)

# anticlust 0.4.0

## Major changes

* `matching()` is a new function for unrestricted or K-partite matching 
to finds groups of similar elements.

* `plot_similarity()`is a new function to plot similarity by cluster
(according to the cluster editing criterion)

* All clustering and anticlustering functions now only take one data 
argument (called `x`) instead of either `features` or `distances`.

* The argument `iv` was removed from `anticlustering()` because it 
does not fit the anticlustering semantic (anticlustering should make
sets «similar» and not dissimilar).

* The random sampling method for anticlustering was removed. 
This implies that the `anticlustering()` function no longer has 
an argument `nrep`.

* The functions `initialize_K()` and `generate_exchange_partners()` were
removed.

* Dropped support for the commercial integer linear programming 
solvers CPLEX and gurobi for exact (anti)cluster editing. If this 
functionality is needed, install version 0.3.0 from Github: 

```
remotes::install_github("m-Py/anticlust", ref = "v0.3.0")
```

* `mean_sd_obj()` no longer computes the discrepancy of 
medians, only in means and standard deviations (as the name also 
suggests).

* In `plot_clusters()`, the arguments `col` and `pch` were removed. 

* In `plot_clusters()`, the argument `clustering` was renamed to `clusters`.

* In `generate_partitions()`, the order of the arguments `N` and 
`K` was switched (the order is now consistent with `n_partitions()`).

* In `balanced_clustering()`, the default `method` was renamed to 
`"centroid"` from `"heuristic"`.

# anticlust 0.3.0

* Release of the package version used in the manuscript 
»Using anticlustering to partition a stimulus pool into equivalent parts«
(Papenberg & Klau, 2019; https://doi.org/10.31234/osf.io/3razc)
