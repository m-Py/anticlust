#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "declarations.h"

/* Exchange Method for weighted cluster editing / clique partitioning
 * 
 * param *data: vector of data points (in R, this is a distance matrix,
 *         the matrix structure must be restored in C)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *frequencies: The number of elements per cluster, i.e., an array
 *         of length *K.
 * param *clusters: An initial assignment of elements to clusters,
 *         array of length *N (MUST BE 1...N)
 * param *order order in which the data is processed, vector of length N.
 * param *mem_error: This is passed with value 0 and only receives the value 1 
 *       if a memory error occurs when executing this function. The caller needs
 *       to test if this value is 1 after execution.
 * 
 * The return value is assigned to the argument `clusters`, via pointer
*/

void wce_heuristic_average(double *data, int *N, int *clusters, int *mem_error) {
        
        const size_t n = (size_t) *N; // number of data points

        // Per data point, store ID, cluster, and category
        struct element POINTS[n]; 
        
        for (size_t i = 0; i < n; i++) {
                POINTS[i].ID = i;
                POINTS[i].cluster = i; // init each element as separate cluster
                POINTS[i].values = malloc(sizeof(double)); // this is just a dummy malloc
                if (POINTS[i].values == NULL) {
                        *mem_error = 1;
                        return;
                }
                POINTS[i].category = 0;
        }
        
        // Some book-keeping variables to track memory error
        int mem_error_cluster_heads = 0;
        int mem_error_cluster_lists = 0;

        /* SET UP CLUSTER STRUCTURE */
        struct node *CLUSTER_HEADS[n];
        mem_error_cluster_heads = initialize_cluster_heads(n, CLUSTER_HEADS);
        
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
                free_cluster_list(CLUSTER_HEADS, n);
                *mem_error = 1;
                return;
        }
        
        // Restore distance matrix
        
        int offsets[n]; // index variable for indexing correct cols in data matrix
        // Allocate memory for distance matrix in C
        double *DISTANCES[n];
        for (size_t i = 0; i < n; i++) {
                DISTANCES[i] = (double*) malloc(sizeof(double) * n);
                if (DISTANCES[i] == NULL) {
                        free_distances(n, DISTANCES, i);
                        free_points(POINTS, n);
                        free_cluster_list(CLUSTER_HEADS, n);
                        *mem_error = 1;
                        return;
                }
        }
        
        // Column offsets (to convert one-dimensional array to Row/Col major)
        for(size_t i = 0; i < n; i++) {
                offsets[i] = i * n;
        }
        
        // Reconstruct the data points as N x N distance matrix
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                        DISTANCES[i][j] = data[offsets[j]++];
                }
        }
        
        
        // Initialize objective
        double OBJ_BY_CLUSTER[n];
        for (size_t i = 0; i < n; i++) {
          OBJ_BY_CLUSTER[i] = 0;
        }

        /* Some variables for bookkeeping during the optimization */
        
        size_t best_cluster;
        double tmp_obj_cl1;
        double tmp_obj_cl2;
        double best_objs[n];
        double tmp_improvement; // encode if an exchange leads to improvement

        // does any improvement occur during exchange method? is used for finding local maximum
        int improvement = 1; 
        
        /* Start main iteration loop for exchange procedure */
        while (improvement == 1) {
                improvement = 0;
                /* 1. Level: Iterate through `n` data points */
                for (size_t i = 0; i < n; i++) {
                        size_t cl1 = PTR_NODES[i]->data->cluster;
                  
                        // Current cluster: Loses distances to element i
                        double sum_dists_i_cl1 = distances_one_element( // can be zero -> good!
                                n, DISTANCES,
                                CLUSTER_HEADS[cl1], i
                        );

                        // Initialize `best` variable for the i'th item
                        double best_improvement = 0;
                        copy_array(n, OBJ_BY_CLUSTER, best_objs);
                        
                        // encode if for element i, it was tested if it should be a singleton cluster
                        int singleton_tried = 0;
                        
                        int exchange_cluster_found = 0;
                        
                        int N_CL1 = number_elements_in_cluster(CLUSTER_HEADS[cl1]);
                        
                        for (size_t u = 0; u < n; u++) {
                                // recode exchange partner index
                                size_t cl2 = u; // OTHER CLUSTER THAT IS TRIED OUT
                                // is it necessary to test the other cluster?
                                if (cl1 == cl2) {
                                        continue;
                                }
                                int N_CL2; 
                                if (CLUSTER_HEADS[cl2]->next == NULL) { // empty cluster
                                        if (singleton_tried == 1) {
                                                continue;
                                        }
                                        N_CL2 = 0;
                                        singleton_tried = 1;
                                } else {
                                  N_CL2 = number_elements_in_cluster(CLUSTER_HEADS[cl2]);
                                }

                                // Initialize `tmp` variables for the exchange partner:
                                tmp_obj_cl1 = OBJ_BY_CLUSTER[cl1]; // These arrays store the sums of distances; averaging is done below
                                tmp_obj_cl2 = OBJ_BY_CLUSTER[cl2];
                                // Update objective
                                tmp_obj_cl1 -= sum_dists_i_cl1;
                                if (N_CL1 > 1) {
                                  tmp_obj_cl1 = tmp_obj_cl1 / (N_CL1 - 1); // average diversity! scale by N
                                }
                                
                                // Other cluster: Gains distances to element i --
                                // compute distances between item i and all items in cluster `cl2`
                                double sum_dists_i_cl2 = distances_one_element(
                                        n, DISTANCES, 
                                        CLUSTER_HEADS[cl2], i
                                );
                                
                                tmp_obj_cl2 += sum_dists_i_cl2;
                                tmp_obj_cl2 = tmp_obj_cl2 / (N_CL2 + 1); // average diversity! scale by N
                                
                                int scale_CL2_before; // prevent division by 0 when averaging diversity
                                if (N_CL2 == 0) {
                                  scale_CL2_before = 1;
                                } else {
                                  scale_CL2_before = N_CL2;
                                }
                                
                                tmp_improvement = (tmp_obj_cl1 + tmp_obj_cl2) -
                                  (OBJ_BY_CLUSTER[cl1] / N_CL1 + OBJ_BY_CLUSTER[cl2] / scale_CL2_before);

                                // Update `best` variables if objective was improved
                                if (tmp_improvement > best_improvement) {
                                        best_objs[cl1] = tmp_obj_cl1;
                                        best_objs[cl2] = tmp_obj_cl2;
                                        best_improvement = tmp_improvement;
                                        best_cluster = cl2;
                                        exchange_cluster_found = 1;
                                }

                        }
                        
                          
                        
                        // Only if objective is improved: Actually *do* the swap
                        if (exchange_cluster_found) {
                                swap_wce(
                                        n, i, cl1, best_cluster,
                                        CLUSTER_HEADS
                                );
                                // Update the "global" variables
                                OBJ_BY_CLUSTER[cl1] = best_objs[cl1];
                                OBJ_BY_CLUSTER[best_cluster] = best_objs[best_cluster];
                                improvement = 1;
                        }
                }
        }


        // Write output
        for (size_t i = 0; i < n; i++) {
                clusters[i] = PTR_NODES[i]->data->cluster;
        }
        
        // in the end, free allocated memory:
        free_points(POINTS, n);
        free_cluster_list(CLUSTER_HEADS, n);
        free_distances(n, DISTANCES, n);
}

int number_elements_in_cluster(struct node *CLUSTER_HEAD) {
  // iterate through current cluster of i
        struct node *tmp = CLUSTER_HEAD;
        int N = 0;
        while (tmp != NULL) {
               N++;
               tmp = tmp->next;
        }
        return N;
}
