
#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

/* Exchange Method for K-Means Anticlustering
 * param *data: vector of data points (in R, this is a data frame / matrix,
 *         but it is passed as a pointer array (vector) to C; here, we use
 *         "clever" indexing on the one-dimensional vector to address the correct 
 *         data)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *M: The number of features (i.e., number of "cols" in *data)
 * param *K: The number of clusters
 * param *frequencies: The number of elements per cluster, i.e., an array
 *         of length *K.
 * param *clusters: An initial assignment of elements to clusters,
 *         array of length *N (has to consists of integers between 0 and (K-1) 
 *         - this has to be guaranteed by the caller)
 * param *partners: A pointer array of length N * k_neighbours, indicating for each 
 *        element which other elements are exchange partners. The first `k_neighbours`
 *        entries are the exchange partners of element 1, the next belong to element
 *        2, and so forth. If an entry is N, the element is skipped (only indices up to N-1 work, 
 *        so N is used to indicate that no exchange partners follow).
 * param *k_neighbours: The number of exchange partners per element.
 * 
 * The return value is assigned to the argument `clusters`, via pointer
 * 
 * 
 * ============================ Some explanations ============================
 * 
 * This is a somewhat different version of kmeans_anticlustering(), which does not
 * use categorical restrictions but can be passed specific exchange partners 
 * for each element (they must be passed to the function and are not generated here).
 * Basically, this function is copy-paste from kmeans_anticlustering() and removes
 * all code using categories and adds code to incorporate custom exchange partners.
 * 
 * ===========================================================================
*/

void fast_kmeans_anticlustering2(double *data, int *N, int *M, int *K, int *frequencies,
                                int *clusters, int *partners, int *k_neighbours, int *mem_error) {

        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of variables per data point
        const size_t k = (size_t) *K; // number of clusters
        const size_t kn = (size_t) *k_neighbours; // total number of exchanges per element
        
        // Some book-keeping variables to track memory error
        int mem_error_points = 0;
        int mem_error_cluster_heads = 0;
        int mem_error_cluster_lists = 0;
        
        // Set up array of data points, fill it, return if memory runs out
        struct element *POINTS;
        POINTS = malloc(n * sizeof(*POINTS)); // free() is included below
        if (POINTS == NULL) { 
                *mem_error = 1;
                return; 
        };

        int USE_CATS = 0;
        int categories = 0;
        
        mem_error_points = fill_data_points( 
                data, n, m, POINTS, clusters, 
                &USE_CATS, &categories
        );

        if (mem_error_points == 1) { // test if all memory was allocated successfully
                free(POINTS);
                *mem_error = 1;
                return;
        }

        /* SET UP CLUSTER STRUCTURE */
        
        // Set up array of pointer-to-cluster-heads
        struct node **CLUSTER_HEADS;
        CLUSTER_HEADS = malloc(k * sizeof(*CLUSTER_HEADS)); // k * pointer to cluster lists
        if (CLUSTER_HEADS == NULL) {
                free_points(POINTS, n);
                free(POINTS);
                *mem_error = 1;
                return; 
        }

        mem_error_cluster_heads = initialize_cluster_heads(k, CLUSTER_HEADS);
        
        if (mem_error_cluster_heads == 1) {
                free(CLUSTER_HEADS);
                free_points(POINTS, n);
                free(POINTS);
                *mem_error = 1;
                return; 
        }

        // Set up array of pointers-to-nodes, return if memory runs out
        struct node **PTR_NODES; // use pointer to pointer as well
        PTR_NODES = malloc(n * sizeof(*PTR_NODES));
        if (PTR_NODES == NULL) {
                free(CLUSTER_HEADS);
                free_points(POINTS, n);
                free(POINTS);
                *mem_error = 1;
                return; 
        }

        mem_error_cluster_lists = fill_cluster_lists(
            n, clusters,
            POINTS, PTR_NODES, CLUSTER_HEADS
        );
        
        if (mem_error_cluster_lists == 1) {
                free_cluster_list(CLUSTER_HEADS, k);
                free(CLUSTER_HEADS);
                free(PTR_NODES);
                free_points(POINTS, n);
                free(POINTS);
                *mem_error = 1;
                return;
        }
        
        // Set up matrix of cluster centers
        double CENTERS[k][m];
        for (size_t i = 0; i < k; i++) {
                compute_center(m, CENTERS[i], CLUSTER_HEADS[i], frequencies[i]);
        }
        
        // Get variance objective of the initial cluster assignment
        double OBJ_BY_CLUSTER[k]; 
        objective_by_cluster(m, k, OBJ_BY_CLUSTER, CENTERS, CLUSTER_HEADS);
        double SUM_OBJECTIVE = array_sum(k, OBJ_BY_CLUSTER);
        
        /* Some variables for bookkeeping during the optimization */
        size_t best_partner;
        double tmp_centers[k][m];
        double best_centers[k][m];
        double tmp_objs[k];
        double best_objs[k];
        double tmp_obj;
        
        /* Start main iteration loop for exchange procedure */
        size_t id_current_exch_partner = 0;
        
        /* 1. Level: Iterate through `n` data points */
        for (size_t i = 0; i < n; i++) {
                size_t cl1 = PTR_NODES[i]->data->cluster;
                
                // Initialize `best` variable for the i'th item
                double best_obj = 0;
                copy_matrix(k, m, CENTERS, best_centers);
                copy_array(k, OBJ_BY_CLUSTER, best_objs);
                
                /* 2. Level: Iterate through the exchange partners */
                for (size_t u = 0; u < kn; u++) {
                        // Get index of current exchange partner
                        size_t j = partners[id_current_exch_partner];
                        // this just counts upwards across all exchange partners, 
                        // ignores matrix-like structure:
                        id_current_exch_partner++;
                        if (j != n) { // no exchange partners any more; coded as index = N
                                size_t cl2 = PTR_NODES[j]->data->cluster;
                                // no swapping attempt if in the same cluster:
                                if (cl1 == cl2) { 
                                  continue;
                                }
                                
                                // Initialize `tmp` variables for the exchange partner:
                                copy_matrix(k, m, CENTERS, tmp_centers);
                                copy_array(k, OBJ_BY_CLUSTER, tmp_objs);
                                
                                update_centers(
                                  k, m, 
                                  tmp_centers, 
                                  PTR_NODES[i],
                                           PTR_NODES[j],
                                                    frequencies
                                );
                                swap(i, j, PTR_NODES);
                                // Update objective
                                tmp_objs[cl1] = cluster_var(m, CLUSTER_HEADS[cl1],
                                                            tmp_centers[cl1]);
                                tmp_objs[cl2] = cluster_var(m, CLUSTER_HEADS[cl2],
                                                            tmp_centers[cl2]);
                                tmp_obj = array_sum(k, tmp_objs);
                                
                                // Update `best` variables if objective was improved
                                if (tmp_obj > best_obj) {
                                  best_obj = tmp_obj;
                                  copy_matrix(k, m, tmp_centers, best_centers);
                                  copy_array(k, tmp_objs, best_objs);
                                  best_partner = j;
                                }
                                
                                // Swap back to test next exchange partner
                                swap(i, j, PTR_NODES);
                        }
                }
                
                // Only if objective is improved: Do the swap
                if (best_obj > SUM_OBJECTIVE) {
                        swap(i, best_partner, PTR_NODES);
                        // Update the "global" variables
                        SUM_OBJECTIVE = best_obj;
                        copy_matrix(k, m, best_centers, CENTERS);
                        copy_array(k, best_objs, OBJ_BY_CLUSTER);
                }
        }

        // Write output
        for (size_t i = 0; i < n; i++) {
                clusters[i] = PTR_NODES[i]->data->cluster;
        }
        
        // in the end, free allocated memory:
        free_points(POINTS, n);
        free_cluster_list(CLUSTER_HEADS, k);
        free(CLUSTER_HEADS);
        free(PTR_NODES);
        free(POINTS);
}
