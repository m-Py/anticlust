
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

Just a vector representing anticluster affiliation

# TODO implementation 

- make function names and argument names consistent
   + data = features (or x?)
   + no item assignment, just anticlustering
   + differentiate functions by objective (distance - variance) rather 
     than by exact - heuristic? not sure though, maybe both is needed
- port functions from minDiff, i.e. repeated random sampling. BUT: use different 
  objectives, i.e. distance and variance objective, abolish self-hacked 
  objective from minDiff
- In addition to repeated random sampling: exchange method / simulated annealing?
  -> both for variance as well as distance criterion
- Tests for all functions. Maybe for ILP function:
    + Test dimensionality of matrices (e.g., constraint matrix)
    + use first version (v 1.0) as benchmark for exact solution; meaning: 
      load all old functions in test files?
- Think about implementing anticlustering algorithms by Valev (1998), even though
  the first I implemented did not seem to really work (see file Stuff.R; this will 
  have to be removed eventually)
- Check: are all dependencies for dplyr and reshape2 removed? Include tests for 
  functions that used these
    + What are the dependencies remaining, just Matrix? How to deal with ILP
      dependencies?
- find out: where did I try out the additional heuristic for distance 
  anticlustering? (the one where I computed max distances and forced these items 
  to be in the same cluster; that was really nice; probably in folder for thesis)

---

- With regard to implementing max-variance anticlustering:

## Try implementing Spaeth's (1987) algorithm

# Uses an "exchange method":

# "That method, applied for (2), improves some (random) initial
# partition by successively moving on trial each object from
# its cluster to all the other ones, and by shifting it, if there is
# any reduction at all, to that one where the first term on the right
# side of (4) decreases most, otherwise taking the next object and
# finally passing through all the objects until no further improvement
# occurs.

