
---

Tests:

- how do all methods deal with missing data?

---

Documentation:

- All exported function fully documented?
- In anticlustering function: use "anticluster editing" as problem 
  description

---

---

`anticlustering` arguments:

- include a verbose argument
    + for complete enumeration and sampling 
      it informs about the progress (x of y iterations; better solution 
      was found etc.)
    + for ILP, print the progress of the solver

--- 

- maybe include a "matching" function
    + this simply wraps "clustering" and "anticlustering" with 
      K = n/2

---

- test if my branch and bound method is useful for clustering 
  (it is not for anticlustering)

--- 

Make interfaces (i.e., arguments) consistent for `balanced_clustering`
and `anticlustering`
