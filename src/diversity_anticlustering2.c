
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "declarations.h"

/* Exchange Method for k-means anticlustering, using pairwise squared Euclidean distances
 * 
 * param *data: vector of data points (N * M features matrix)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *M: The number of features (i.e., number of "columns" in *data)
 * param *clusters: An initial assignment of elements to clusters,
 *         array of length *N (MUST BE 1...N)
 * param *mem_error: This is passed with value 0 and only receives the value 1 
 *       if a memory error occurs when executing this function. The caller needs
 *       to test if this value is 1 after execution.
 * 
 * The return value is assigned to the argument `clusters`, via pointer
*/

void diversity_kmeans_anticlustering(double *data, int *N, int *M, int *K, int *clusters, int *mem_error) {
        
        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of features
        const size_t k = (size_t) *K; // number of clusters

        // Per data point, store ID, cluster, and category
        struct element POINTS[n]; 
        // Create offset variable to correctly read out data points
        int m_ptr[m];
        for (size_t i = 0; i < m; i++) {
              m_ptr[i] = i * n;
        }
        
        for (size_t i = 0; i < n; i++) {
                POINTS[i].ID = i;
                POINTS[i].cluster = clusters[i]; // init each element as separate cluster
                POINTS[i].values = malloc(sizeof(double) * m); // this is just a dummy malloc
                // Fill data into `element`:
                if (POINTS[i].values == NULL) {
                        *mem_error = 1;
                        return;
                }
                for (size_t j = 0; j < m; j++) {
                        POINTS[i].values[j] = data[m_ptr[j]++];
                }
                POINTS[i].category = 0;
        }
        
        // Some book-keeping variables to track memory error
        int mem_error_cluster_heads = 0;
        int mem_error_cluster_lists = 0;

        /* SET UP CLUSTER STRUCTURE */
        struct node *CLUSTER_HEADS[k];
        mem_error_cluster_heads = initialize_cluster_heads(k, CLUSTER_HEADS);
        
        if (mem_error_cluster_heads == 1) {
                free_points(POINTS, n);
                *mem_error = 1;
                return; 
        }

        // Set up array of pointers-to-nodes, return if memory runs out
        struct node *PTR_NODES[n];
        mem_error_cluster_lists = fill_cluster_lists(
                n, clusters,
                POINTS, PTR_NODES, CLUSTER_HEADS
        );
        if (mem_error_cluster_lists == 1) {
                free_points(POINTS, n);
                free_cluster_list(CLUSTER_HEADS, k);
                *mem_error = 1;
                return;
        }
        
        // Initialize objective
        double OBJ_BY_CLUSTER[n];
        for (size_t i = 0; i < n; i++) {
                OBJ_BY_CLUSTER[i] = 0;
        }

        /* Some variables for bookkeeping during the optimization */
        
        size_t best_partner;
        size_t best_cluster;
        double tmp_obj_cl1;
        double tmp_obj_cl2;
        double best_objs[k];
        double tmp_improvement; // encode if an exchange leads to improvement
        
        // does any improvement occur during exchange method? is used for finding local maximum
        int improvement = 1; 
        
        /* Start main iteration loop for exchange procedure */
        while (improvement == 1) {
                improvement = 0;
                /* 1. Level: Iterate through `n` data points */
                for (size_t i = 0; i < n; i++) {
                        size_t cl1 = PTR_NODES[i]->data->cluster;
                        double* features_i = PTR_NODES[i]->data->values;
                  
                        // Current cluster: Loses distances to element i
                        double sum_dists_i_cl1 = squared_distances_one_element(
                                n, m, features_i, CLUSTER_HEADS[cl1]
                        );

                        // Initialize `best` variable for the i'th item
                        double best_improvement = 0;
                        copy_array(k, OBJ_BY_CLUSTER, best_objs);
                        int exchange_cluster_found = 0;
                        
                        for (size_t u = 0; u < n; u++) {

                                size_t cl2 = PTR_NODES[u]->data->cluster;
                                double* features_j = PTR_NODES[u]->data->values;
                                // is it necessary to test the other cluster?
                                if (cl1 == cl2) {
                                        continue;
                                }
                                
                                
                               
                                // Initialize `tmp` variables for the exchange partner:
                                tmp_obj_cl1 = OBJ_BY_CLUSTER[cl1];
                                tmp_obj_cl2 = OBJ_BY_CLUSTER[cl2];

                                // CL1 loses distances to i, gains distances to j
                                // CL2 loses distances to j, gains distances to i
                                double sum_dists_j_cl2 = squared_distances_one_element(
                                        n, m, features_j, CLUSTER_HEADS[cl2]
                                );
                                double sum_dists_i_cl2 = squared_distances_one_element(
                                      n, m, features_i, CLUSTER_HEADS[cl2]
                                ) - euclidean_squared(features_i, features_j, m); // they cannot be in the same cluster!
                                double sum_dists_j_cl1 = squared_distances_one_element(
                                  n, m, features_j, CLUSTER_HEADS[cl1]
                                ) - euclidean_squared(features_i, features_j, m);
                                
                                
                                // Update objective
                                tmp_improvement = -sum_dists_i_cl1 - sum_dists_j_cl2 + sum_dists_i_cl2 + sum_dists_j_cl1;
                                tmp_obj_cl1 += -sum_dists_i_cl1 + sum_dists_j_cl1; 
                                tmp_obj_cl2 += -sum_dists_j_cl2 + sum_dists_i_cl2; 
                                
                                // Update `best` variables if objective was improved
                                if (tmp_improvement > best_improvement) {
                                        best_objs[cl1] = tmp_obj_cl1;
                                        best_objs[cl2] = tmp_obj_cl2;
                                        best_improvement = tmp_improvement;
                                        best_cluster = cl2;
                                        best_partner = u;
                                        exchange_cluster_found = 1;
                                }
                                
                        }
                        // Only if objective is improved: Actually *do* the swap
                        if (exchange_cluster_found) {
                                swap(i, best_partner, PTR_NODES);
                                // Update the "global" variables
                                OBJ_BY_CLUSTER[cl1] = best_objs[cl1];
                                OBJ_BY_CLUSTER[best_cluster] = best_objs[best_cluster];
                                improvement = 0; // do not use local maximum search, just exchange method 
                        }
                }
        }


        // Write output
        for (size_t i = 0; i < n; i++) {
                clusters[i] = PTR_NODES[i]->data->cluster;
        }
        
        // in the end, free allocated memory:
        free_points(POINTS, n);
        free_cluster_list(CLUSTER_HEADS, k);
}

// compute sum of distances of one node to other nodes in the same cluster
double squared_distances_one_element(size_t n, size_t m, double* features, struct node *start_node) {
        struct node *tmp = start_node->next;
        double sum = 0;
        while (tmp != NULL) {
                // self distance will be 0 so it is not treated as special case
                sum += euclidean_squared(features, tmp->data->values, m); 
                tmp = tmp->next;
        }
        return sum;
}
