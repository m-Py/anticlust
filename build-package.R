pckDIR <- "/home/martin/git/anticlust/"
myRlib <- "/home/martin/R/x86_64-pc-linux-gnu-library/3.6/"
#myRlib <- "/home/martin/R/i686-pc-linux-gnu-library/3.3/"

library("roxyPackage")

pckg.dscrptn <- data.frame(
  Package="anticlust",
  Type = "Package",
  Title = "Subset Partitioning via Anticlustering",
  AuthorsR = "c(person(given='Martin', family='Papenberg',
       email='martin.papenberg@hhu.de', role=c('aut', 'cre', 'cph')), 
person(given='Meik', family='Michalke', role=c('ctb'), comment='centroid based 
clustering algorithm')
     )",
  Description="Anticlustering partitions a pool of elements into subsets 
(i.e., anticlusters) in such a way that the subsets are as similar as 
possible. This is accomplished by maximizing instead of minimizing a 
clustering objective, such as the intra-cluster variance or the sum of 
pairwise distances within clusters. This way, anticlustering maximizes the 
data heterogeneity within all clusters. The package includes exact and 
heuristic approaches for anticlustering. The exact approach requires that a 
linear programming solver is installed. The package also offers methods for 
classical clustering problems under size constraints, i.e., creating balanced 
clusters of equal size.",
  License="MIT + file LICENSE",
  Encoding="UTF-8",
  LazyData = "true",
  URL="https://github.com/m-Py/anticlust",
  stringsAsFactors = FALSE,
  Suggests = "gurobi, Rglpk, Rcplex, testthat, knitr, rmarkdown, dplyr",
  Imports = "Matrix, RANN",
  Depends = "R (>= 3.4.0)",
  BugReports = "https://github.com/m-Py/anticlust/issues",
  VignetteBuilder = "knitr, rmarkdown"
)

roxy.package(
   pck.source.dir = pckDIR,
   pck.version = "0.3.0-1",
   R.libs = myRlib,
   repo.root = "~/R/repo/anticlust",
   pck.description = pckg.dscrptn,
   actions = c("roxy", "package", "cite", "html", "doc", "cleanRd", "buildVignettes") # "check"
) 
