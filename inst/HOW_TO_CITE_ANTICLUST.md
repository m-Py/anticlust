# Citing the anticlust package

If you use any anticlustering methods from the `anticlust` package in your research, it would be courteous if you cite the following paper that introduces the package:

- Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition data sets into equivalent parts. *Psychological Methods, 26*(2), 161--174. https://doi.org/10.1037/met0000301

The functions that implement anticlustering methods are currently:

- `anticlustering()`
- `kplus_anticlustering()`
- `bicriterion_anticlustering()`
- `fast_anticlustering()`

If you use the bicriterion anticlustering heuristic by Brusco et al. (2020), you should also cite their paper: 

- Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining diversity and dispersion criteria for anticlustering: A bicriterion approach. *British Journal of Mathematical and Statistical Psychology, 73*(3), 375–396. https://doi.org/10.1111/bmsp.12186

You are using the algorithm by Brusco et al. (2020) if you use the function `anticlustering()` with the argument `method = "brusco"` or if you are using the function `bicriterion_anticlustering()`. These functions implement the combined bicriterion algorithm MBPI + BILS. 

If you use k-plus anticlustering, you should cite the following reference:

- Papenberg, M. (2023). K-plus Anticlustering: An Improved k-means Criterion for Maximizing Between-Group Similarity. *British Journal of Mathematical and Statistical Psychology*. Advance online publication. https://doi.org/10.1111/bmsp.12315

You are using k-plus anticlustering if you use the function `anticlustering()` with the argument `objective = "kplus"` or if you are using the function `kplus_anticlustering()`.

## References for anticlustering algorithms

The default anticlustering algorithm `method = "exchange"` was described in Papenberg and Klau (2021). However, this algorithm is actually just a shorter version of the local maximum method (i.e., `method = "local-maximum"`), which was used in Papenberg (2023) and corresponds to the algorithm "LCW" by Weitz and Lakshminarayanan (1996, 1998). By default, the anticlustering algorithms (including LCW) create equal-sized groups, and Papenberg (2023) discussed how to adapt the LCW method for different sized groups. However, this extension is pretty straight forward and has probably been done elsewhere. 

Papenberg and Klau (2021) presented an adaptation of the LCW method to incorporate categorical constraints in the anticlustering process (i.e,, ensure that a nominal variable is as evenly distributed as possible between groups). This extension is used when specifying the `categories` argument in `anticlustering()`. The function `categorical_sampling()` implements the stratified split that initially balances the categories across clusters, and is internally called by `anticlustering()`. Papenberg and Klau (2021) also discuss using a reduced number of exchange partners when the LCW algorithm runs; this happens when using the function `fast_anticlustering()` and specifying the argument `k_neighbours`.

The optimal anticlustering method based on integer linear programming (i.e., `method = "ilp"`) was presented in Papenberg and Klau (2021). This ILP is an extension of the model by Grötschel and Wakabayashi (1989).

## References for anticlustering objectives

A classical reference discussing the **diversity** objective is Feo and Khellaf (1990). Brusco et al. (2020) were the first to discuss the diversity in the context of Psychological research. Papenberg and Klau (2021) and Papenberg (2023) conducted a comparative evaluation of the k-means and diversity objective (using the Euclidean distance to compute the diversity). See Brusco et al. (2020) and Papenberg and Klau (2021) for additional references on the diversity objective (e.g., Gallego et al., 2013). Many other papers have discussed algorithms for maximizing the diversity (i.e., solving the Maximum Diverse Grouping Problem, MDGP), even very recently (e.g. Yang et al., 2022). Note that the reversal of the diversity, which is the sum of intra-cluster (dis)similarities, is a classical clustering objective that has already been discussed in the literature much earlier (e.g., Zahn, 1964).

Späth (1986) discussed the reversal of the **k-means** criterion (i.e., the **"variance"**) and coined the term anticlustering for this purpose. Note that Valev (1983, 1998) also independently and earlier introduced the term anticlustering. Reviews of the enormousness body of literature on k-means *clustering* are for example found in Jain (2010) and Steinley (2006).

In the context of anticlustering problems, the **dispersion** objective has for example been discussed by Fernández et al. (2013) and Brusco et al. (2020). However note that the reversal of the dispersion is also a classical clustering objective and much earlier work exists, for example in the context of hierarchical clustering (see references in Brusco et al. 2020).

**K-plus** anticlustering was introduced by Papenberg (2023). The k-plus criterion extends the k-means objective, but it no longer has a (reversed) clustering interpretation. Thus, no earlier references exist that discuss the k-plus criterion in the context of cluster analysis. 

This list on literature on anticlustering objectives is necessarily not exhaustive. However, if you feel that important references are missing, please tell me <a href="mailto:martin.papenberg@hhu.de">via email</a>.

## Cluster analysis

Note that the following `anticlust` functions do not implement "anti"-clustering methods, but instead clustering / matching methods (i.e., the intra-group similarity is high and the between-group similarity is low):

- `wce()`
- `balanced_clustering()`
- `matching()`

With the exception of `balanced_clustering(..., method = "ilp")`, it may not be entirely appropriate to cite Papenberg and Klau (2021) when using these functions, because that paper is concerned with "anti"-clustering methods. However, it would be nice if you just cite the R package `anticlust` itself, e.g. via:

Papenberg, M. (2019). Anticlust: Subset partitioning via anticlustering (Version 0.6.4) [Computer software]. https://cran.r-project.org/package=anticlust

For `wce()`, the "scientific" reference would be Grötschel and Wakabayashi (1989), because it directly implements their model. For other possible references, see `?wce`. 

The functions `balanced_clustering(..., method = "centroid)` (i.e., the default method) and `matching()` implement a centroid (anti)clustering algorithm developed by [Meik Michalke](https://www.psychologie.hhu.de/arbeitsgruppen/diagnostik-und-differentielle-psychologie/arbeitsgruppe). The details of this algorithm have not (yet?) been published in a paper, and its functionality is currently only described in the `anticlust` documentation, see `?balanced_clustering`. For now, you should probably stick with citing the `anticlust` package itself when using the functions `balanced_clustering()` or `matching()`. That, or you find an appropriate way of citing Meik's work using a named reference, e.g. by referring to the commit where the centroid method was contributed to `anticlust`:  [https://github.com/m-Py/anticlust/commit/e5e7ab2d452aa](https://github.com/m-Py/anticlust/commit/e5e7ab2d452aa). When using the function `balanced_clustering()` with the argument `method = "ilp"`, you can safely cite Papenberg and Klau (2021) because this method has been described in the section "Preclustering". 

## References

Brusco, M. J., Cradit, J. D., & Steinley, D. (2020). Combining diversity and dispersion criteria for anticlustering: A bicriterion approach. *British Journal of Mathematical and Statistical Psychology, 73*(3), 375–396. https://doi.org/10.1111/bmsp.12186

Feo, T. A., & Khellaf, M. (1990). A class of bounded approximation algorithms for graph partitioning. *Networks, 20*, 181--195.

Fernández, E., Kalcsics, J., & Nickel, S. (2013). The maximum dispersion problem. *Omega, 41*, 721--730. https://doi.org/10.1016/j.omega.2012.09.005

Gallego, M., Laguna, M., Marti, R., & Duarte, A. (2013). Tabu search with strategic oscillation for the
maximally diverse grouping problem. *Journal of the Operational Research Society, 64*, 724--
734. https://doi.org/10.1057/jors.2012.66

Grötschel, M., & Wakabayashi, Y. (1989). A cutting plane algorithm for a clustering problem. *Mathematical Programming, 45*, 59-96.

Jain, A. K. (2010). Data clustering: 50 years beyond k-means. *Pattern Recognition Letters, 31*, 651--666.

Papenberg, M., & Klau, G. W. (2021). Using anticlustering to partition data sets into equivalent parts. *Psychological Methods, 26*(2), 161--174. https://doi.org/10.1037/met0000301

Papenberg, M. (2023). K-plus Anticlustering: An Improved k-means Criterion for Maximizing Between-Group Similarity. *British Journal of Mathematical and Statistical Psychology*. Advance online publication. https://doi.org/10.1111/bmsp.12315

Steinley, D. (2006). K-means clustering: A half-century synthesis. *British Journal of Mathematical and Statistical Psychology, 59*, 1--34.

Valev, V. (1983). Set partition principles. In J. Kozesnik (Ed.), *Transactions of the ninth Prague conference on information theory, statistical decision functions, and random processes* (Prague, 1982) (pp. 251--256). Prague, Czech Republic: Springer Netherlands.

Valev, V. (1998). Set partition principles revisited. In A. Amin, D. Dori, P. Pudil, & H. Freeman (Eds.), *Advances in pattern recognition, SSPR/SPR, Lecture notes in computer science* (pp. 875– 881). Heidelberg,
Germany: Springer.

Weitz, R. R., & Lakshminarayanan, S. (1996). On a heuristic for the final exam scheduling problem. *Journal of the Operational Research Society, 47*(4), 599--600.

Weitz, R. R., & Lakshminarayanan, S. (1998). An empirical comparison of heuristic methods for creating maximally diverse groups. *Journal of the Operational Research Society, 49*(6), 635--646. https://doi.org/10.1057/palgrave.jors.2600510

Yang, X., Cai, Z., Jin, T., Tang, Z., & Gao, S. (2022). A three-phase search approach with dynamic population size for solving the maximally diverse grouping problem. *European Journal of Operational Research, 302*(3), 925--953.

Zahn, C. T. (1964). Approximating symmetric relations by equivalence relations. *Journal of the Society for Industrial and Applied Mathematics, 12*, 840--847.

