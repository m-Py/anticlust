# anticlust 0.2.5

2019-05-27

The new version of anticlust now enables parallelization of the random 
sampling method, improving run time.

- The `anticlustering` function now has an additional argument 
  (`parallelization`) that can be used to activate parallel computation 
  when using the heuristic method
- For now, the default value of `parallelization` is `FALSE`
- Another argument was added (`seed`) to make the random sampling method 
  reproducible
    + Just using `set.seed` prior to the computation does not make 
    the function call reproducible when `parallelization` is `TRUE` 
    because each core has its own random seed

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
