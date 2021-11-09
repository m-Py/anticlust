# anticlust 0.6.0

* `anticlust` now includes the bicriterion algorithm for simultaneously 
maximizing diversity and dispersion, proposed by Brusco et al.
(<doi:10.1111/bmsp.12186>) and implemented by Martin Breuer (for details see [his bachelor thesis](https://www.cs.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Algorithmische_Bioinformatik/Bachelor-_Masterarbeiten/2516084_ba_ifo_AbschlArbeit_klau_mapap102_mabre121_20200820_2320.pdf))
  * It can be called from the main function `anticlustering()` by setting 
  `method = "brusco"`; in this case only either dispersion or diversity is
  maximized
  * `bicriterion_anticlustering()` -- newly exported in this version -- can be used 
  for a more fine grained usage of the Brusco et al. algorithm, fully using
  its main functionality to optimize both dispersion as well as diversity

# anticlust 0.5.7

- Just an update to the documentation: Updating all references to the Papenberg & Klau paper after its "actual" publication in Psychological Methods.

# anticlust 0.5.6-1

- Minor bug fix in `plot_clusters()` (via [87f585798](https://github.com/m-Py/anticlust/commit/87f5857986d16b338827a213ea0e95dda8e86eef))

# anticlust 0.5.6

### User-visible changes

- `plot_clusters()` now uses the default color palette to highlight 
the different clusters 
- `plot_clusters()` now uses different `pch` symbols when the number of 
clusters is low (K < 8)

### Internal changes

* `anticlustering()` and `categorical_sampling()` now better balance categorical 
variables when the output groups require different sizes (i.e., if the group
sizes do not have any common denominator)

* Some additional input validations for more useful error messages when arguments
in `anticlustering()` are not correctly specified

# anticlust 0.5.5

## New feature

* `anticlustering()` has a new argument `standardize` to standardize the data
input before the optimization starts. This is useful to give all variables the 
same weight in the anticlustering process, irregardless of the scaling of the
variables. Especially useful for `objective = "kplus"` to ensure that both 
minimizing differences with regard to means and variance is equally important.

## Bug fix

* Fixes a memory leak in the C code base, via [2c4fe6d](https://github.com/m-Py/anticlust/commit/2c4fe6d4c272b10717158e639a76103ef574cc41)

# anticlust 0.5.4-1

- Internal change: `anticlustering()` with `objective = "dispersion"` now 
<<<<<<< HEAD
implements the local updating procedure [proposed by Martin Breuer](https://www.cs.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Algorithmische_Bioinformatik/Bachelor-_Masterarbeiten/2516084_ba_ifo_AbschlArbeit_klau_mapap102_mabre121_20200820_2320.pdf). 
=======
implements the local updating procedure [proposed by Martin Breuer](https://www.cs.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Algorithmische_Bioinformatik/Bachelor-_Masterarbeiten/2516084_ba_ifo_AbschlArbeit_klau_mapap102_mabre121_20200820_2320.pdf).
>>>>>>> devel
This leads to a considerable speedup when maximizing the dispersion, enabling the fast 
processing of large data sets.

# anticlust 0.5.4

## User-visible changes

- `anticlustering()` now has native support for the maximizing the dispersion
  objective, setting `objective = "dispersion"`. The dispersion is the minimum 
  distance between any two elements within the same cluster, see 
  `?dispersion_objective`.

## Internal changes 

- The exchange optimization algorithm for anticlustering has been reimplemented 
  in C, leading to a substantial boost in performance when using one of the 
  supported objectives "diversity", "variance", "dispersion", or "kplus".
  (Optimizing user-defined objective functions still has to be done in plain 
  R and therefore has not been sped up.)

# anticlust 0.5.3

* `kplus_objective()` is a new function to compute the value of the k-plus
  criterion given a clustering. See `?kplus_objective` for details.

* In `anticlustering()` and `categorical_sampling()`, the argument 
  `K` can now be used to specify the size of the groups, not just the 
  number of groups. This way, it is easy to request groups of different
  size. See the help pages `?anticlustering` and `?categorical_sampling`
  for examples.

# anticlust 0.5.2-1 / 0.5.2-2

* Fixed two minor bugs that prevented the correct transformation of class `dist` to 
  class `matrix` when using the repeated exchange (or "local-maximum") method, 
  see [c42e136](https://github.com/m-Py/anticlust/commit/c42e1367ec371dc054a5dd51916b45e1424d6274) 
  and [e6fdae5](https://github.com/m-Py/anticlust/commit/e6fdae50965150781d1f4621844f24c63167364a).

# anticlust 0.5.2

## User-visible changes

* In `anticlustering()`, there is a new option for the argument `method`:
  "local-maximum". When using `method = "local-maximum"`, the exchange method is 
  repeated until an local maximum is reached. That means after the exchange 
  process has been conducted for each data point, the algorithm restarts with 
  the first element and proceeds to conduct exchanges until the objective cannot 
  be improved. This procedure is more in line with classical neighbourhood search
  that only terminates when a local optimum is reached.

* In `anticlustering()`, there is now a new argument `repetitions`. It can be used
  to specify the number of times the exchange procedure (either `method = 
  "exchange"` or `method = "local-maximum"`) is called. `anticlustering()`
  returns the best partitioning found across all repetitions.
  
* `anticlustering()` now implements a new objective function, extending the classical
  k-means criterion, given by `objective = "kplus"`. Using `objective = 
  "kplus"` will minimize differences with regard to both means and standard deviations 
  of the input variables, whereas k-means only focuses on the means. Details on this 
  objective will follow.  

# anticlust 0.5.1

- Fixes a bug in `anticlustering()`, that led to an incorrect 
computation of cluster centers with option `objective = "variance"` 
for unequal cluster sizes, see
[2ef6547](https://github.com/m-Py/anticlust/commit/2ef65475d9516d93f7a1950f3e3af30a561e52da)

# anticlust 0.5.0

## User-visible changes

### Major

* A new exported function: `categorical_sampling()`. Categorical 
sampling can be used to obtain a stratified split of a data set. Using 
this function is like calling `anticlustering()` with argument 
`categories`, but no clustering objective is maximized. The categories 
are just evenly split between samples, which is very fast (in contrast 
to the exchange optimization that may take some time for large data 
sets). Apart from the categorical restriction that balances the 
frequency of categories between samples, the split is random.

* The function `distance_objective()` was renamed into 
`diversity_objective()` because there are several clustering objectives 
based on pairwise distances, e.g. see the new function 
`dispersion_objective()`.

* `dispersion_objective()` is a new function to compute the dispersion 
of a given clustering, i.e., the minimum distance between two elements 
within the same group. Maximizing the dispersion is an anticlustering 
task, see the help page of `dispersion_objective()` for an example.

### Minor

* Several changes to the documentation, in particular now highlighting 
the publication of the paper "Using Anticlustering to Partition 
Data Sets Into Equivalent Parts" (https://doi.org/10.1037/met0000301) 
describing the algorithms and criteria used in the package `anticlust` 

* In `anticlustering()`, anticluster editing is now by default requested 
using `objective = "diversity"` (but `objective = "distance"` is still 
supported and leads to the same behaviour). This change was done because
there are several anticlustering objectives based on pairwise distances.

* `anticlustering()` can no longer use an argument `K` of length > 1 
with `preclustering = TRUE` because this resulted in undocumented 
behaviour (this is a good change because it does not make sense to 
specify an initial assignment of elements to groups via `K` and at the 
same time request that preclustering handles the initial assignment)

* When using a custom objective function, the order of the required 
arguments is now reversed: The data comes first, the clustering second.

* Because the order of arguments in custom objective functions was 
reversed, the function `mean_sd_obj()` now has reversed arguments as 
well.

* The package vignettes are no longer distributed with the package 
itself because rendering R Markdown resulted in an error with the 
development version of R. This may change again in the future when R 
Markdown no longer throws an error with R devel. The vignette is 
currently available via the package website 
(https://m-py.github.io/anticlust/).

## Internal changes

* Improved running speed of generating constraints in integer linear
programming variant of (anti)clustering, via 
[0a870240f8](https://github.com/m-Py/anticlust/commit/0a870240f8264f0e74f4cbf0b20d789cfa0d6469)

# anticlust 0.4.1

## User-visible changes

* In `anticlustering()`, preclustering and categorical constraints can 
now be used at the same time. In this case, exchange partners are 
clustered within the same category, using a call to `matching()` passing 
`categories` to argument `match_within`.

* In `anticlustering()`, it is now possible to use `preclustering = 
TRUE` for unbalanced data size (e.g., if N = 9 and K = 2).

* In `matching()`, it is now possible to prevent sorting the output by 
similarity using a new argument `sort_output`. Its default is `TRUE`, 
setting it to `FALSE` prevents sorting. This prevents some extra 
computation that is necessary to determine similarity for each cluster.

## Minor 

* Some changes to documentation

* There is now a package website at https://m-py.github.io/anticlust/

* Additional error handling

## Internal

* Improvements to implementation of k-means anticlustering (i.e., in 
`anticlustering()` with `objective == "variance"` or in 
`fast_anticlustering()`)
  * on each exchange iteration, only recomputes distances from clusters 
whose elements have been swapped (improves run time relevant for larger 
K).
  * Previously, there were only as many exchange partners per element as 
members in the least frequent category if argument `categories` was 
passed). This was not documented behavior and is undesirable. Now, all 
members from a category may serve as exchange partners, even if the 
categories have different size.

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
