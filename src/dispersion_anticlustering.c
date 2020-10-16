
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "declarations.h"

/* Exchange Method for Anticlustering Based on a Distance matrix
 * 
 * param *data: vector of data points (in R, this is a distance matrix,
 *         the matrix structure must be restored in C)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *K: The number of clusters
 * param *frequencies: The number of elements per cluster, i.e., an array
 *         of length *K.
 * param *clusters: An initial assignment of elements to clusters,
 *         array of length *N (has to consists of integers between 0 and (K-1) 
 *         - this has to be guaranteed by the caller)
 * int *USE_CATS A boolean value (i.e., 1/0) indicating whether categorical 
 *         constraints
 * param *C: The number of categories
 * param *CAT_frequencies: The number of elements per category, i.e., an array
 *         of length *C.
 * param *categories: An assignment of elements to categories,
 *         array of length *N (has to consists of integers between 0 and (C-1) 
 *         - this has to be guaranteed by the caller)
 * param *mem_error: This is passed with value 0 and only receives the value 1 
 *       if a memory error occurs when executing this function. The caller needs
 *       to test if this value is 1 after execution.
 * 
 * The return value is assigned to the argument `clusters`, via pointer
*/

void dispersion_anticlustering(double *data, int *N, int *K, int *clusters, 
                              int *USE_CATS, int *C, int *CAT_frequencies,
                              int *categories, int *mem_error) {
        
        const size_t n = (size_t) *N; // number of data points
        const size_t k = (size_t) *K; // number of clusters
        
        // Per data point, store ID, cluster, and category
        struct element POINTS[n]; 
        
        for (size_t i = 0; i < n; i++) {
                POINTS[i].ID = i;
                POINTS[i].cluster = clusters[i];
                POINTS[i].values = malloc(sizeof(double)); // this is just a dummy malloc
                if (POINTS[i].values == NULL) {
                        *mem_error = 1;
                        return;
                }
                if (*USE_CATS) {
                        POINTS[i].category = categories[i];
                } else {
                        POINTS[i].category = 0;
                }
        }
        
        // Some book-keeping variables to track memory error
        int mem_error_categories = 0;
        int mem_error_cluster_heads = 0;
        int mem_error_cluster_lists = 0;
        
        // Deal with categorical restrictions
        size_t c = number_of_categories(USE_CATS, C);
        *CAT_frequencies = get_cat_frequencies(USE_CATS, CAT_frequencies, n);
        
        size_t *CATEGORY_HEADS[c];
        mem_error_categories = get_indices_by_category(
                n, c, CATEGORY_HEADS, USE_CATS, categories, CAT_frequencies, POINTS
        );
        if (mem_error_categories == 1) {
                free_points(n, POINTS, n);
                *mem_error = 1;
                return; 
        }
        
        /* SET UP CLUSTER STRUCTURE */
        struct node *CLUSTER_HEADS[k];
        mem_error_cluster_heads = initialize_cluster_heads(k, CLUSTER_HEADS);
        
        if (mem_error_cluster_heads == 1) {
                free_points(n, POINTS, n);
                free_category_indices(c, CATEGORY_HEADS, c);
                *mem_error = 1;
                return; 
        }

        // Set up array of pointers-to-nodes, return if memory runs out
        struct node *PTR_NODES[n];
        mem_error_cluster_lists = fill_cluster_lists(
                n, k, clusters, 
                POINTS, PTR_NODES, CLUSTER_HEADS
        );
        
        if (mem_error_cluster_lists == 1) {
                free_points(n, POINTS, n);
                free_category_indices(c, CATEGORY_HEADS, c);
                free_cluster_list(k, CLUSTER_HEADS, k);
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
                        free_points(n, POINTS, n);
                        free_category_indices(c, CATEGORY_HEADS, c);
                        free_cluster_list(k, CLUSTER_HEADS, k);
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
        double OBJECTIVE = dispersion_objective(n, k, DISTANCES, CLUSTER_HEADS);

        /* Some variables for bookkeeping during the optimization */
        
        size_t best_partner;
        double tmp_obj;
        
        /* Start main iteration loop for exchange procedure */
        
        /* 1. Level: Iterate through `n` data points */
        for (size_t i = 0; i < n; i++) {
                size_t cl1 = PTR_NODES[i]->data->cluster;
                
                // Initialize `best` variable for the i'th item
                double best_obj = 0;
                
                /* 2. Level: Iterate through the exchange partners */
                size_t category_i = PTR_NODES[i]->data->category;
                // `category_i = 0` if `USE_CATS == FALSE`.
                size_t n_partners = CAT_frequencies[category_i];
                // `CAT_frequencies[0] == n` if `USE_CATS == FALSE`
                for (size_t u = 0; u < n_partners; u++) {
                        // recode exchange partner index
                        size_t j = CATEGORY_HEADS[category_i][u];
                        size_t cl2 = PTR_NODES[j]->data->cluster;
                        // no swapping attempt if in the same cluster:
                        if (cl1 == cl2) { 
                                continue;
                        }
                        
                        // Using local updating of dispersion objective
                        int i_has_dispersion = has_node_dispersion(
                                n, DISTANCES, OBJECTIVE, 
                                CLUSTER_HEADS[PTR_NODES[i]->data->cluster], 
                                PTR_NODES[i], 0
                        );
                        int j_has_dispersion = has_node_dispersion(
                                n, DISTANCES, OBJECTIVE, 
                                CLUSTER_HEADS[PTR_NODES[j]->data->cluster], 
                                PTR_NODES[j], 0
                        );
                        
                        int before = i_has_dispersion || j_has_dispersion;
                        
                        // Now swap the elements in the cluster list
                        swap(n, i, j, PTR_NODES);
                        
                        i_has_dispersion = has_node_dispersion(
                                n, DISTANCES, OBJECTIVE, 
                                CLUSTER_HEADS[PTR_NODES[i]->data->cluster], 
                                PTR_NODES[i], 1
                        );
                        j_has_dispersion = has_node_dispersion(
                                n, DISTANCES, OBJECTIVE, 
                                CLUSTER_HEADS[PTR_NODES[j]->data->cluster], 
                                PTR_NODES[j], 1
                        );
                        
                        int after = i_has_dispersion || j_has_dispersion;
                        
                        if (after) {
                                tmp_obj = best_obj - 1; // cannot improve objective
                        } else if (!before) {
                                tmp_obj = best_obj - 1; // cannot improve objective
                        } else if (before) {
                                tmp_obj = dispersion_objective(n, k, DISTANCES, CLUSTER_HEADS);
                                // todo: more speed up
                        }

                        // Update `best` variables if objective was improved
                        if (tmp_obj > best_obj) {
                                best_obj = tmp_obj;
                                best_partner = j;
                        }
                        // Swap back to test next exchange partner
                        swap(n, i, j, PTR_NODES);
                }
                
                // Only if objective is improved: Do the swap
                if (best_obj > OBJECTIVE) {
                        swap(n, i, best_partner, PTR_NODES);
                        // Update the "global" variables
                        OBJECTIVE = best_obj;
                }
        }

        // Write output
        for (size_t i = 0; i < n; i++) {
                clusters[i] = PTR_NODES[i]->data->cluster;
        }
        
        // in the end, free allocated memory:
        free_points(n, POINTS, n);
        free_category_indices(c, CATEGORY_HEADS, c);
        free_cluster_list(k, CLUSTER_HEADS, k);
        free_distances(n, DISTANCES, n);
}

// Check if an element has a connection to another element in the same cluster, where 
// their pairwise distance is the current dispersion. For implementing Martin B.'s clever speed ups.
// Has parameter `after` to check if the minimum within-cluster distance is lower after 
// a swap.
int has_node_dispersion(size_t n, double *distances[n], double dispersion, struct node *HEAD, struct node *node, int after) {
        struct node *current = HEAD->next;
        int id_node = node->data->ID;
        while (current != NULL) {
                int tmp_id = current->data->ID;
                current = current->next;
                if (tmp_id == id_node) {
                        continue; 
                }
                double current_distance = distances[id_node][tmp_id];
                if (after) {
                        if (current_distance <= dispersion) {
                                return 1;
                        }
                } else {
                        if (current_distance == dispersion) {
                                return 1;
                        }
                }
                
        }
        return 0;
}


// Compute sum of distances by cluster
double dispersion_objective(size_t n, size_t k, double *distances[n], 
                          struct node *HEADS[k]) {
        double min = INFINITY;
        for (size_t i = 0; i < k; i++) {
                double distance = minimun_distance_cluster(n, distances, HEADS[i]);
                if (distance < min) {
                        min = distance;
                }
        }
        return min;
}

// Compute minimum distance within a cluster
double minimun_distance_cluster(size_t n, double *distances[n], struct node *HEAD) {
        double min = INFINITY;
        struct node *current = HEAD->next;
        while (current != NULL) {
                double distance = minimin_distance_element(n, distances, current, current->data->ID);
                if (distance < min) {
                        min = distance;
                }
                current = current->next;
        }
        return min;
}

// compute minimum distance of one element to other elements in the same cluster
double minimin_distance_element(size_t n, double *distances[n], 
                                struct node *start_node, size_t ID) {
        struct node *tmp = start_node->next;
        double min = INFINITY;
        while (tmp != NULL) {
                double distance = distances[ID][tmp->data->ID];
                if (distance < min) {
                        min = distance;
                }
                tmp = tmp->next;
        }
        return min;
}
