
#include "declarations.h"
#include <stdlib.h> 

/* Extracted method that fill data points into array of struct `element`
 * param `double *data` pointer to original data array of length n
 * param `size_t m`: Number of data points
 * param `size_t m`: Number of variables per data point
 * param `struct element points[n]`: Array to be filled with data points
 * param `int *clusters`: Cluster affiliation of the n data points
 * 
 * return: `0` if all data points were successfully stored; `1` if not.
 * 
 */
int fill_data_points(double *data, size_t n, size_t m, struct element *POINTS, 
                     int *clusters, int *USE_CATS, int *categories) {
        // Create offset variable to correctly read out data points
        int m_ptr[m];
        for (size_t i = 0; i < m; i++) {
                m_ptr[i] = i * n;
        }
        
        // Size of a data vector per element:
        size_t data_size = m * sizeof(POINTS[0].values[0]);

        for (size_t i = 0; i < n; i++) {
                POINTS[i].cluster = clusters[i];
                if (*USE_CATS) {
                        POINTS[i].category = categories[i];
                } else {
                        POINTS[i].category = 0;
                }
                POINTS[i].ID = i;
                
                POINTS[i].values = malloc(data_size);
                if (POINTS[i].values == NULL) {
                        free_points(POINTS, i);
                        return 1;
                } 
                
                // Fill data into `element`:
                for (size_t j = 0; j < m; j++) {
                        POINTS[i].values[j] = data[m_ptr[j]++];
                }

        }
        return 0;
}

/* After creation, initialize the HEAD of each cluster list
 * 
 * We access the cluster list using an array of length `k` pointing
 * to the `k` HEADs of the clusters. The first element pointed to 
 * in the list (i.e., the HEAD) is "empty". It just points to the 
 * first real element that is in that cluster.
 * 
 * param `size_t k`: The number of clusters
 * param `struct node *PTR_CLUSTER_HEADS[k]`: The array of pointers to 
 *     cluster HEADS.
 * 
 *  * return: `0` if the cluster list could be initialized successfully, `1` 
 *      if not (in that case, there was no memory that could be allocated).
 * 
 */

int initialize_cluster_heads(size_t k, struct node **HEADS) {
        for (size_t i = 0; i < k; i++) {
                HEADS[i] = (struct node*) malloc(sizeof(struct node));
                if (HEADS[i] == NULL) {
                        free_cluster_list(HEADS, i);
                        return 1;
                }
                HEADS[i]->next = NULL;
        }
        return 0;
}

/* After initialization, fill the cluster lists with data
 * This function does two things at the same time:
 * (a) add each data point as a node to a cluster list,
 * (b) store the pointer to each node in the array `PTR_NODES`
 */
int fill_cluster_lists(size_t n, int *clusters,
                       struct element *POINTS, struct node **PTR_NODES,
                       struct node **PTR_CLUSTER_HEADS) {
        for (size_t i = 0; i < n; i++) {
                struct node *CLUSTER_HEAD = PTR_CLUSTER_HEADS[clusters[i]];
                PTR_NODES[i] = append_to_cluster(CLUSTER_HEAD, &POINTS[i]);
                if (PTR_NODES[i] == NULL) { // failed to allocate memory
                        return 1;
                }
        }
        return 0;
}

/* Append data point to linked list
 * 
 * param `struct node *HEAD` Pointer to HEAD of cluster list
 * param `struct element *data`: Pointer to the data point 
 *     that is appended to the list
 * 
 * return: Pointer to the `node` of the element that was appended to the
 *     cluster list
 * 
 */

struct node* append_to_cluster(struct node *HEAD, struct element *data) {
        struct node *tmp = HEAD->next; // may be NULL if list is empty
        HEAD->next = (struct node*) malloc(sizeof(struct node));
        if (HEAD->next == NULL) {
                return NULL; // failed to allocate memory
        }
        // New element is right next to HEAD element:
        HEAD->next->data = data;
        HEAD->next->next = tmp; // tmp was next to HEAD before
        return HEAD->next; 
}
