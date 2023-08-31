
#include <stdlib.h> 
#include "declarations.h"
#include <stdio.h>

/* Free memory in the cluster lists
* param `size_t k`: The maximum number of clusters
* param `struct node *PTR_CLUSTER_HEADS[k]`: The array of pointers to 
*     cluster HEADS
* param size_t i: The actual number of clusters; usually the same 
*     as `k`, but may be different if setting up the cluster list 
*     fails during setup. (the same logic holds for the below free 
*     functions that target different data structures)
*/
void free_cluster_list(struct node **PTR_CLUSTER_HEADS, size_t i) {
        struct node *ptr;
        struct node *prev; // using temp pointer for freeing
        for (size_t j = 0; j < i; j++) {
                ptr = PTR_CLUSTER_HEADS[j];
                while (ptr->next != NULL)
                {  
                        prev = ptr;
                        ptr = ptr->next;
                        free(prev);
                }
                free(ptr);
        }
}

/* Free index array for categories */
void free_category_indices(size_t **CATEGORY_HEADS, size_t i) {
    for (size_t j = 0; j < i; j++) {
        free(CATEGORY_HEADS[j]);
    }
}

void free_distances(size_t n, double *DISTANCES[n], size_t i) {
        for (size_t j = 0; j < i; j++) {
            free(DISTANCES[j]);
        }
}

/* Free memory in the data points
 * param `size_t n`: length of array `POINTS`
 * param `struct element POINTS[n]`: Array containing data points
 * param `size_t i`: The index up to which data points are freed
 */
void free_points(struct element *POINTS, size_t i) {
        for (size_t j = 0; j < i; j++) {
                free(POINTS[j].values);
        }
}
