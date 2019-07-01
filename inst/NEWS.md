
# anticlust 0.2.7

2019-07-01

A big update:

- A new algorithm for anticlustering is available: the exchange method.
  It is based on investigating the "neighborhood" of a given partition and
  swapping the anticlusters of two elements. The neighborhood of an 
  element is defined as all possible swaps that can be done with elements 
  that are currently assigned to a different cluster. The swap is made that 
  improves the objective function the most across all possible swaps.
  This method outperforms the random sampling heuristic for k-means 
  anticlustering but it has quadratic run time; when setting 
  `preclustering = TRUE`, the number of swaps is reduced because only
  very similar elements are legal swap partners, making the exchange
  method applicable for larger N. The exchange method may incorporate
  both categorical and preclustering constraints, which is not possible
  for the random sampling approach. For anticluster editing, 
  the random sampling approach is better than the exchange method.
  As there are now two heuristic methods (random sampling and exchange) 
  the argument `method` of the function `anticlustering()`
  now has the following three possible values: "sampling", "exchange", 
  "ilp". In earlier versions, the two options were "heuristic" and "ilp"; 
  hence this change possibly breaks earlier code.
- A new function `generate_partitions` can be used to generate all
  partitions, making it possible to solve anticlustering via complete
  enumeration. In particular, it is now possible---for small
  problem instances---to solve k-means anticlustering optimally, which
  cannot be done with integer linear programming.  The help file 
  (`?generate_partitions`) contains example code illustrating how to do 
  this.

Minor changes:

- The "variance" criterion can now be computed when there are missing 
  values in the input data.
- In `plot_clusters`, it is now possible to adjust the size of the 
  cluster centroid using the new argument `cex_centroid`
  

# anticlust 0.2.6

2019-06-19

Minor update: `plot_clusters` now has an additional argument 
`illustrate_variance`. If this argument is set to `TRUE`, a cluster
solution is illustrated with the k-means variance criterion.

# anticlust 0.2.5

2019-05-27

The new version of anticlust now enables parallelization of the random 
sampling method, improving run time.

- The `anticlustering` function now has an additional argument 
  (`parallelize`) that can be used to activate parallel computation 
  when using the heuristic method
- For now, the default value of `parallelize` is `FALSE`
- Another argument was added (`seed`) to make the random sampling method 
  reproducible
    + Just using `set.seed()` prior to the computation does not make 
    a function call reproducible when `parallelize` is `TRUE` 
    because each core has its own random seed
    + The `seed` argument is optional

An example data set is now included with the package, courteously 
provided by Marie Lusia Schaper and Ute Bayen. For details, see 
`?schaper2019`.

# anticlust 0.2.4

2019-04-26

The new version of anticlust includes support for constraints induced 
by grouping variables. 

- The `anticlustering` function now has an optional argument (`categories`)
  that can be used to induce categorical constraints
- `categories` can be a vector if there is is one grouping variable or 
  a matrix/data if there is more than one grouping variable
- Currently, `categories` can only be used with the random sampling 
  method (`method = "heuristic"`)
- `categories` overrides the value of `preclustering`; it is not 
  possible to use categorical and preclustering constraints at the same
  time

In `anticlustering`, the default value of `preclustering` is now FALSE.
