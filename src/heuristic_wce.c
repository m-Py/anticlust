
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

void wce_heuristic(double *data, int *N, int *clusters, int *mem_error) {
        
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
        double SUM_OBJECTIVE = 0;

        /* Some variables for bookkeeping during the optimization */
        
        size_t best_cluster;
        double tmp_objs[n];
        double best_objs[n];
        double tmp_obj;
        
        /* Start main iteration loop for exchange procedure */
        int improvement = 1;

        while (improvement == 1) {
                improvement = 0;
                /* 1. Level: Iterate through `n` data points */
                for (size_t i = 0; i < n; i++) {
                        size_t cl1 = PTR_NODES[i]->data->cluster;
                        
                        // Initialize `best` variable for the i'th item
                        double best_obj = 0;
                        copy_array(n, OBJ_BY_CLUSTER, best_objs);
                        
                        for (size_t u = 0; u < n; u++) {
                                // recode exchange partner index
                                size_t cl2 = u; // OTHER CLUSTER THAT IS TRIED OUT
                                // Initialize `tmp` variable for the exchange partner:
                                copy_array(n, OBJ_BY_CLUSTER, tmp_objs);
                                // Update objective
                                // Current cluster: Loses distances to element i
                                tmp_objs[cl1] -= distances_one_element( // can be zero -> good!
                                        n, DISTANCES,
                                        CLUSTER_HEADS[cl1], i
                                );
                                
                                // Swap!
                                swap_wce(
                                        n, i, cl1, cl2, 
                                        PTR_NODES, CLUSTER_HEADS
                                );
                                
                                // Other cluster: Gains distances to element 1
                                tmp_objs[cl2] += distances_one_element(
                                        n, DISTANCES, 
                                        CLUSTER_HEADS[cl2], i
                                );

                                tmp_obj = array_sum(n, tmp_objs);
                                
                                // Update `best` variables if objective was improved
                                if (tmp_obj > best_obj) {
                                        best_obj = tmp_obj;
                                        copy_array(n, tmp_objs, best_objs);
                                        best_cluster = cl2;
                                }
                                // Swap back to test next exchange partner
                                swap_wce(
                                        n, i, cl2, cl1,
                                        PTR_NODES, CLUSTER_HEADS
                                );
                        }
                        
                        // Only if objective is improved: Do the swap
                        if (best_obj > SUM_OBJECTIVE) {
                                swap_wce(
                                        n, i, cl1, best_cluster,
                                        PTR_NODES, CLUSTER_HEADS
                                );
                                // Update the "global" variables
                                SUM_OBJECTIVE = best_obj;
                                copy_array(n, best_objs, OBJ_BY_CLUSTER);
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

// swap procedure for swapping with different cluster label, not with element in other cluster
// i = current item
// cl2 = cluster with which it is swapped
int swap_wce(size_t n, size_t i, size_t cl1, size_t cl2, struct node *PTR_NODES[n], struct node *CLUSTER_HEADS[n]) {
        // iterate through current cluster of i
        struct node *tmp = CLUSTER_HEADS[cl1];
        while (tmp != NULL) {
                if (tmp->next && tmp->next->data->ID == i) {
                        struct node *COPY = tmp->next;
                        tmp->next = tmp->next->next; // REMOVAL
                        // Add to new cluster, next to HEAD; store previous next of HEAD
                        struct node *COPY2 = CLUSTER_HEADS[cl2]->next;
                        CLUSTER_HEADS[cl2]->next = COPY; // ADDED
                        CLUSTER_HEADS[cl2]->next->data->cluster = cl2; // UPDATED CLUSTER
                        CLUSTER_HEADS[cl2]->next->next = COPY2; // LINKED TO PREVIOUS LIST
                        return 0; // job's done
                }
                tmp = tmp->next;
        }
        return 0; // should not come here
}
