
# How to form the API

## Functions

Do I have one function for group assignment or do I offer several
functions? In the first case, the algorithm would be selected through a
parameter; in the latter case, each algorithm gets one function.

- Pro one function
    + for the user, it is easy to see how the package is
      used maybe, 
    + maybe, it is possible to set good defaults that do not have to be
      changed often
- Pro several functions
    + several different parameters have to be adjusted and not just the
      one that selects algorithm (e.g., for repeated random assignment,
      it is necessary to specify the number of assignments; this does
      not make any sense for exact or heuristic approaches)

**Initial approach**

- Make several functions that have common parameter usage (i.e., make
  function usage as consistent as possible for different methods);
  describe these functions clearly in documentation.
- Later: check if it even makes sense to abstract a single function from
  all methods
- Maybe: mix the two approaches: some algorithms may be better suited to
  be stored in one function, others may need a separate function

**What functions do I export**

1. Exact ILP approach 
2. Heuristic ILP (based on preclustering)
3. More levels of "heuristicism" 
    + Exact preclustering, repeated random assignment
    + Heuristic preclustering, repeated random assignment
4. Repeated random assignment (port from minDiff)
    + Employ a different default objective: the one that is also used in
      the other functions
    + It should remain possible to pass your own functions as the
      objective
    + Rename function create_groups to a more descriptive name that
      should describe the employed algorithm, maybe: "enumerate_groups"
5. anticlustering [to be implemented]
6. Algorithm for distribution-balanced stratified cross-validation [to
   be implemented]
7. Cluster editing algorithms
    + Weighted cluster editing
    + p-cluster editing
    + min/max cluster editing
    + equal-sized cluster editing

**What will be the "default" function?**

What function will be suggested to the user at first?

## Return value

Do I return a `data.frame` where that contains the group assignment as
an additional column (this is currently implemented in most functions
and also in the package minDiff) or only a vector? Maybe a simple vector
is the way to go that can be easily implemented consistently across
functions.

# TODO implementation 

- tidy up function `item_assign_ilp`
- implement the constraint ILP matrix as sparse matrix
- do not create the ILP twice for heuristic assignment 
- make consistent return value for the different assignment methods
- port functions from minDiff
- Function `obj_value` should take two input parameters: one for the
  assignment, and another one that represents the item distances
