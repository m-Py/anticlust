
#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

/* Exchange Method for Anticlustering
 * param *data: vector of data points (in R, this is a data frame,
 *         the matrix structure must be restored in C)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *M: The number of features (i.e., number of "cols" in *data)
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
 * 
 * 
 * ============================ Some explanations ============================
 * 
 * Throughout this method, some data structures are defined that are used
 * throughout. The function works by first initializing the required data 
 * structures on which an exchange optimization algorithm is concuted.
 * 
 * These are the data structures that are "global" throughout the function:
 * 
 * 1. `struct element POINTS[n]` - An array copy of the `n` data points in the 
 *    same order as the input. Stores the `m` values per data point and the 
 *    cluster assignment of the data point (initially passed via `*clusters`)
 * 2. `struct node *CLUSTER_HEADS[k]` - An array of length `k` where each
 *    entry points to the head of a list that represents a cluster. Each cluster
 *    is implemented as a linked list where each node points to an `element` in 
 *    `struct element POINTS[n]`. Thus, there are `k` cluster lists that are 
 *    used during the algorithm to compute the variance by cluster. During
 *    the exchange method, elements are swapped between cluster lists 
 *    (implemented through the function `swap()`). 
 *    The function `initialize_cluster_heads()` sets up the pointer array; the 
 *    function `fill_cluster_lists()` fills the data points as nodes into the 
 *    lists.
 * 3. `struct node *PTR_NODES[n]` - Array of pointers to nodes, used for 
 *    iterating during the exchange method. Points to elements in the cluster 
 *    lists but can be used to iterate through the data in the original order
 *    of the data (across cluster lists, this order is lost).
 * 4. `size_t *CATEGORY_HEADS[c]` - Array of pointers to indices. 
 *    The i'th element of `C_HEADS` contains an array of all 
 *    indices of the elements that have the i'th category.
 * 5. `double CENTERS[k][m]` - A matrix of cluster centers.
 * 6. `double OBJECTIVE_BY_CLUSTER[k]` - The variance objective by cluster.
 * 7. `double SUM_VAR_OBJECTIVE` - The total variance objective (i.e., sum 
 *    of all entries in `OBJECTIVE_BY_CLUSTER`).
 *    
 * Regarding the inclusion of categorical constraints: If the argument 
 * `USE_CATS` is `TRUE` (i.e, `1`), the exchange method restricts the 
 * exchange partners to elements of the same category. What elements are part 
 * of the same category is specified via the input argument `categories`. 
 * Generally, the arguments `C`, `CAT_frequencies`, `categories` have the same 
 * semantic as `K`, `frequencies` and `clusters`, respectively. However, they 
 * represent a fixed category and not a variable cluster affiliation (that is 
 * changed in this function to maximize the k-means criterion). If `USE_CATS` is
 * `FALSE` (i.e, `0`), the arguments `C`, `CAT_frequencies`, `categories` are 
 * not used.
 * 
 * To balance a categorical variable that is represented by `categories` across
 * clusters, it is necessary that these are already balanced when calling this
 * function. It will not obtain an initial balanced partitioning itself - 
 * the caller is responsible.
 * 
 * ===========================================================================
*/

void kmeans_anticlustering(double *data, int *N, int *M, int *K, int *frequencies,
        int *clusters, int *USE_CATS, int *C, int *CAT_frequencies,
        int *categories, int *mem_error) {
        
        /* 
        * - Free strategy: each function cleans its own mess
        *   + free in low-level function when possible (i.e., when an allocation error 
        *     occurs within this function)
        *   + Free in high-level function when previously, memory was allocated successfully,
        *     but then an allocation error occured
        *   + TODO use a safe-free function instead of just `free()`
        */ 
    
        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of variables per data point
        const size_t k = (size_t) *K; // number of clusters
        
        // Some book-keeping variables to track memory error
        int mem_error_points = 0;
        int mem_error_categories = 0;
        int mem_error_cluster_heads = 0;
        int mem_error_cluster_lists = 0;
        
        // Set up array of data points, fill it, return if memory runs out
        struct element POINTS[n];
        mem_error_points = fill_data_points(
                data, n, m, POINTS, clusters, 
                USE_CATS, categories
        );
        
        if (mem_error_points == 1) {
                *mem_error = 1;
                return;
        }
        
        /* CATEGORICAL RESTRICTIONS 
         * Write the pointer of arrays `CATEGORY_HEADS`
        */
        
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
        
        // Set up array of pointer-to-cluster-heads
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
        
        /* 1. Level: Iterate through `n` data points */
        for (size_t i = 0; i < n; i++) {
                size_t cl1 = PTR_NODES[i]->data->cluster;
                
                // Initialize `best` variable for the i'th item
                double best_obj = 0;
                copy_matrix(k, m, CENTERS, best_centers);
                copy_array(k, OBJ_BY_CLUSTER, best_objs);
                
                /* 2. Level: Iterate through the exchange partners */
                size_t category_i = PTR_NODES[i]->data->category;
                size_t n_partners = CAT_frequencies[category_i];
                for (size_t u = 0; u < n_partners; u++) {
                        // Get index of current exchange partner
                        size_t j = CATEGORY_HEADS[category_i][u];
                        
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
                        swap(n, i, j, PTR_NODES);
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
                        swap(n, i, j, PTR_NODES);
                }
                
                // Only if objective is improved: Do the swap
                if (best_obj > SUM_OBJECTIVE) {
                        swap(n, i, best_partner, PTR_NODES);
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
        free_points(n, POINTS, n);
        free_category_indices(c, CATEGORY_HEADS, c);
        free_cluster_list(k, CLUSTER_HEADS, k);
}

/* 
 * Perform a swap between two elements.
 * 
 * param `size_t n`: Number of data points
 * param `size_t i`: Index of first element to be swapped
 * param `size_t j`: Index of second element to be swapped
 * param `struct node *PTR_NODES[n]`: pointer to nodes
 * 
 */

void swap(size_t n, size_t i, size_t j, struct node *PTR_NODES[n]) {
        
        struct node *one = PTR_NODES[i];
        struct node *two = PTR_NODES[j];
        
        // Get cluster indices
        size_t cl1 = one->data->cluster;
        size_t cl2 = two->data->cluster;
        
        // Update pointers to elements
        size_t ID1 = one->data->ID;
        size_t ID2 = two->data->ID;
        PTR_NODES[ID1] = two;
        PTR_NODES[ID2] = one;
        
        // Update the cluster affiliation
        one->data->cluster = cl2;
        two->data->cluster = cl1;
        
        // Update nodes in cluster lists 
        struct element *tmp = one->data;
        one->data = two->data;
        two->data = tmp;
        
}
