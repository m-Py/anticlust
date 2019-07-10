
# anticlust 0.2.9

2019-07-09

- A bug was fixed that led to an incorrect computation of the objective 
  function for anticluster editing when employing the exchange method (see 
  [243ca64](https://github.com/m-Py/anticlust/commit/243ca642be787e8c59ece4dbbb1b567fdac05656)). 
  Tests show that the exchange method now outperforms random sampling
  for anticluster editing (as well as for k-means anticlustering). 
  Therefore, the exchange method is now the default method 
  (see [b101073](https://github.com/m-Py/anticlust/commit/b101073602906b6b9bbf00c76943668f43407e0e)).
- The fast exchange method is now used when optimizing the variance 
  criterion in a call to `anticlustering()`. This improves run time by 
  a large margin for this important application. See [2f47fea](https://github.com/m-Py/anticlust/commit/2f47feaf05aee1d53b60bf78bb7c02994a4659c9).
- Two changes with regard to the functionality of arguments in `anticlustering()`
    + It is possible that the argument `objective` now takes as input a 
      function. The passed function has to take two arguments,
      the first being a cluster assignment vector (such as returned by 
      `anticlustering()`), the second being the data the objective is 
      computed on (e.g. an N x M matrix where rows are elements and 
      columns are features). Larger return values must indicate a 
      better objective as the objective is maximized with the existing 
      methods (exchange method and random sampling). This functionality 
      makes it possible for users to implement
      their own operationalization of set similarity.
    + It is possible that the argument `preclustering` now takes as input
      a preclustering vector and not only `TRUE` or `FALSE` (in the former 
      case, the preclustering vector has been computed within the 
      `anticlustering()` function). This allows for more flexibility 
      in combining preclustering and anticlustering methods. For 
      example, it is now possible to conduct optimal preclustering using 
      integer linear programming with the function `balanced_clustering()`,
      and then use a heuristic anticlustering method that incorporates this 
      preclustering.
    + Both of these changes have not yet been added to the function 
      documentation as they require some more testing.
- The `fast_anticlustering()` function has been documented more
  thoroughly and part of the `anticlustering()` docs have been reworked
  (now advocating the exchange method as the preferable option).

# anticlust 0.2.8

2019-07-05

A new exported function is available: `fast_anticlustering()`. As the 
name suggests, it is optimized for speed and particularly useful for large 
data sets (many thousand elements). It uses the k-means variance 
objective because computing all pairwise distances for the cluster 
editing objective becomes computationally infeasible for large data
sets. Additionally, it employs a speed-optimized exchange method. 
The number of exchange partners can be adjusted using the argument 
`k_neighbours`. Fewer exchange partners make it possible to apply the
`fast_anticlustering()` function to very large data sets. The default
value for `k_neighbours` is `Inf`, meaning that in the default case, 
each element is swapped with all other elements.

# anticlust 0.2.7

2019-07-01

A big update:

- A new algorithm for anticlustering is available: the exchange method.
  Building on an initial random assignment, elements are swapped between 
  anticlusters in such a way that each swap improves anticluster similarity by 
  the largest amount that is possible (cf. Späth, 1986).
  This procedure is repeated for each element; because each 
  possible swap is investigated for each element, the total number of 
  exchanges grows quadratically with input size, rendering the exchange 
  method unsuitable for large N. Setting `preclustering = TRUE` will 
  limit the legal exchange partners to very similar elements, resulting 
  in improved run time while preserving a rather good solution.
  The exchange method outperforms the random sampling heuristic for k-means 
  anticlustering. The exchange method may incorporate
  both categorical and preclustering constraints, which is not possible
  for the random sampling approach. As there are now two heuristic 
  methods (random sampling and exchange) 
  the argument `method` of the function `anticlustering()`
  now has the following three possible values: "sampling", "exchange", 
  "ilp". In earlier versions, the two options were "heuristic" and "ilp"; 
  this change does not break earlier code because using 
  `method = "heuristic"` will still refer to the random sampling method.
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
