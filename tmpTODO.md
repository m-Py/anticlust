
---

Tests:

- how do all methods deal with missing data?

---

Documentation:

- All exported function fully documented?

---


Implementation:

- function `anticlustering` should also take distance object as input
- complete enumeration as alternative to ILP for method = "exact" if 
  no solver is installed 

---

`anticlustering` arguments:

- include a verbose argument
    + for complete enumeration and sampling 
      it informs about the progress (x of y iterations; better solution 
      was found etc.)
    + for ILP, print the progress of the solver

- Rework method argument:
    + "exact" is used for the exact computation only; 
      it checks whether an ILP solver is available, if YES, the ILP 
      is used, otherwise complete enumeration is used
    + Add "ilp" and "enumeration" methods that specificall ask for 
      integer linear programming and complete enumeration. Preclustering 
      can be activated for "ilp"
      
--- 

- maybe include a "matching" function
    + this simply wraps "clustering" and "anticlustering" with 
      K = n/2

---

- test if the branch and bound method is useful for clustering 
