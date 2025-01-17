# anticlust 0.8.9

- `categories_to_binary()` no longer uses dummy coding with a reference category, but instead codes each levels of a categorical variable with a separate variable (thanks to [Gunnar Klau](https://github.com/guwek) for spotting potential problems with dummy coding).
- `anticlustering()` has new `method = "2PML"`, which is an improved heuristic when using must-link constraints

# anticlust 0.8.8

## User visible changes

- `three_phase_search_anticlustering()` implements the three phase search algorithm by [Yang et al.](https://doi.org/10.1016/j.ejor.2022.02.003), contributed by Hannah Hengelbrock ([@HanneyAI](https://github.com/HanneyAI))

# anticlust 0.8.7

## User visible changes

- `anticlustering()` now has an argument `must_link`, which can be used force elements into the same cluster
- It is now possible that a cluster (e.g., in `anticlustering()`) only has 1 member (this threw an error before)

## Bug fixes

- `diversity_objective()` is now computed correctly when a cluster only has one member (fixed via [8403fab1461b2cda8](https://github.com/m-Py/anticlust/commit/8403fab1461b2cda8dda768d8e80c5ec92552e4a))
- Fixed a memory leak in `anticlustering(..., objective = "diversity")` thanks to [@HanneyAI](https://github.com/HanneyAI) (via [24c244faf8b2c0774071](https://github.com/m-Py/anticlust/commit/24c244faf8b2c07740710cac957a753f96f545fc))

# anticlust 0.8.6

## New features 

- `bicriterion_anticlustering()` has new arguments: `dispersion_distances`, `average_diversity`, `init_partitions`, `return`.
- `anticlustering()` now has new `objective = "average-diversity"`
- In `anticlustering()`, `method = "brusco"` now works for `objective = "variance"` and `"objective = kplus"`
- Added `lpSolve` solver as backend for the optimal methods, and it is now the default solver
- `optimal_anticlustering()` and `optimal_dispersion()` now have an additional argument `time_limit`
- `anticlustering()` now has an argument `cannot_link`, which can be used to forbid pairs of elements from being assigned to the same cluster. When used, this solves the same (NP hard) graph coloring problem as `optimal_dispersion()`. Unlike the other optimal methods, it uses the Symphony solver with priority, when it is found (otherwise the `lpSolve`)

## Bug fixes

- Bug fix in `optimal_dispersion()`: Output element `$edges` no longer includes edges that were investigated in the last iteration of the algorithm (and which are not relevant for finding the optimal dispersion)

## Internal changes / Other 

- Speed improvements for `anticlustering(..., objective = "diversity")` when using `method = "local-maximum"` and `repetitions` (the restart algorithm is now entirely implemented in C and does not call `method = "exchange"` repeatedly from R)
- `anticlust` now depends on package `lpSolve`

# anticlust 0.8.5

## New features

- `optimal_anticlustering()` is a new exported function that gathers all currently (and in the future) implemented optimal algorithms for anticlustering
- `balanced_clustering()` now has an argument `solver`, which can be used to specify the ILP solver when using `method = "ilp"`

## Internal changes

- The default selection of ILP solvers in `anticlustering()`, `balanced_clustering()` and `optimal_dispersion()` was changed due to a reoccurring CRAN issue: If both the Rglpk and the Rsymphony packages are available, the GLPK will now be prioritized. This is because the SYMPHONY solver sometimes crashes on Macs (or at least on one CRAN test station). The `optimal_anticlustering()`, `optimal_dispersion()`, and `balanced_clustering()` functions have an argument `solver` that can be used to circumvent this default behaviour.
- `anticlust` now uses [`tinytest`](https://cran.r-project.org/package=tinytest) instead of [`testthat`](https://cran.r-project.org/package=testthat) for unit tests.

# anticlust 0.8.3

## Documentation

- Some minor updates to documentation and vignettes
- Updating all references to the k-plus anticlustering paper after its "actual" publication:

Papenberg, M. (2024). K-plus Anticlustering: An Improved k-means Criterion for Maximizing Between-Group Similarity. *British Journal of Mathematical and Statistical Psychology, 77*(1), 80--102. https://doi.org/10.1111/bmsp.12315

## Internal changes

- `fast_anticlustering()` received [another internal change](https://github.com/m-Py/anticlust/commit/8d3c85dc0076) to improve the speed of the re-computation of the objective during the optimization. In particular, updating the objective is now done by only inspecting the two clusters between which an exchange actually took place, instead of re-computing a sum across all clusters.

# anticlust 0.8.2

## Internal changes

- (Regression) `anticlustering(..., objective = "variance")` uses pre 0.8.0 implementation to fix some CRAN issues

# anticlust 0.8.0

## New features

- `fast_anticlustering()` now has an additional argument `exchange_partners`, which can be used to pass custom exchange partners instead of using the default nearest neighbour search.
-  `generate_exchange_partners()` is a new exported function that can be used to address the new argument `exchange_partners` in `fast_anticlustering()`.

## Internal changes 

- `anticlustering()` received internal changes to ensure that it [no longer crashes the computer for about N > 250000 elements](https://github.com/m-Py/anticlust/issues/50).
- `fast_anticlustering()` has been re-implemented in C, which is much faster than the previous R implementation.
- `fast_anticlustering()` now uses an alternative computation of the k-means objective, which reduces run time by an order of magnitude as compared to before.

## Documentation 

- Expanded documentation of `fast_anticlustering()`.
- The vignette "Speeding up anticlustering" has been rewritten to reflect that `fast_anticlustering()` is now again the best choice for processing (very) large data sets.

# anticlust 0.7.0

- An exact ILP method is now available for maximizing the dispersion, contributed by Max Diekhoff.
  * `optimal_dispersion()` is a new exported function implementing the method
  * `anticlustering()` makes it available when using `method = "ilp"` and `objective = "dispersion"`
-  `kplus_moment_variables()` is a new exported function that generates k-plus variables from a data set
  * Offers some additional flexibility as compared to calling `kplus_anticlustering()`, which generates these variables internally (e.g., use k-plus augmentation on some variables but not all -- such as binary variables)
- `categories_to_binary()` is a new exported function that converts one or several categorical variables into a binary representation
  * Can be used to include categorical variables as part of the optimization criterion in k-means / k-plus anticlustering, see new vignette "Using categorical variables with anticlustering"
- 3 new vignettes have been added to the `anticlust` documentation
- Fixed a bug in `kplus_anticlustering()` that did not correctly implement `preclustering = TRUE`
- It is now possible to use the SYMPHONY solver as backend for the optimal ILP methods.

# anticlust 0.6.4-1 / 0.6.4-2 / 0.6.4-3

- [Implements](https://github.com/m-Py/anticlust/commit/435958a0577) 
[some](https://github.com/m-Py/anticlust/commit/9bb8275cad) 
[fixes](https://github.com/m-Py/anticlust/commit/33cee784cb392a) in the 
internal function `gdc_set()` that finds the greatest common denominator in a 
set of numbers. The fixes prevent `categorical_sampling()` (which is also called by `anticlustering()` when using the `categories` argument) from potentially running 
into an infinite loop when combining uneven group sizes via `K` with a 
`categories` argument.

# anticlust 0.6.4

- `kplus_anticlustering()` now has an argument `T` instead of `moments`, where `T` denotes the number of distribution moments considered during k-plus anticlustering (`moments` was an integer vector specifying each individual moment that should be considered)
  * Explanation: Lower order moments should be skipped in favour of higher order moments, so the new interface makes more sense.

# anticlust 0.6.3

**Major changes**

* This release adds a new exported function and removes two others (I very much doubt anyone used those, though -- see below -- if your code is affected, please email me).
  - `kplus_anticlustering()` is a new exported function: A new interface function to k-plus anticlustering, implementing the k-plus method as described in "K-plus Anticlustering: An Improved K-means Criterion for Maximizing Between-Group Similarity" (Papenberg, 2023; https://doi.org/10.1111/bmsp.12315). Using `anticlustering(x, K, objective = "kplus")` is still supported and remains unchanged. The new function `kplus_anticlustering()`, however, offers more functionality and nuance with regard to optimizing the k-plus objective family.
  - The function `kplus_objective()` was removed. 
  - The function `mean_sd_obj()` was removed. 
  
Explanations for the rather drastic changes, i.e., removing instead of deprecating functions (that very likely do not affect anyone):

- Given the advanced theoretical background for k-plus anticlustering, the function  `kplus_objective()`  no longer makes any sense. Given that the k-plus objective is a family of objectives, keeping the function that computes one special case is more harmful to keep it than to just remove it now. As the k-plus objective basically re-uses the k-means criterion, maintaining a function such as `kplus_objective()` was questionable to begin with.

- Since there is the k-plus anticlustering method now, I did not want to keep the "hacky" way to optimize similarity with regard to means and standard deviations, i.e., using the `mean_sd_obj()` function as `objective` in anticlustering. Please use the k-plus method to optimize similarity with regard to means and standard deviations (you can even extend to skewness, kurtosis, and other higher order moments; see the new `kplus_anticlustering()` function).

**Minor changes**

- Finally added Marie Luisa Schaper as contributor for contributing her data set
- Some work on documentation

# anticlust 0.6.2

- Some work on docs and examples

# anticlust 0.6.1

- Minor bug fix in C code base via [c1a5604f](https://github.com/m-Py/anticlust/commit/c1a5604f36a08cb963428af46d98ece616c73568)

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

- Just an update to the documentation: Updating all references to the Papenberg & Klau paper after its "actual" publication in Psychological Methods:

Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition data sets into equivalent parts. *Psychological Methods, 26*(2), 161–174. https://doi.org/10.1037/met0000301

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
implements the local updating procedure [proposed by Martin Breuer](https://www.cs.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Algorithmische_Bioinformatik/Bachelor-_Masterarbeiten/2516084_ba_ifo_AbschlArbeit_klau_mapap102_mabre121_20200820_2320.pdf). 
This leads to a considerable speedup when maximizing the dispersion, enabling the fast processing of large data sets.

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
