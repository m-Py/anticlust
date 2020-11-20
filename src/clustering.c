#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

/* 
 Function for balanced clustering, implementing centroid algorithm by Meik Michalke
 
 Expects that the order in which neighbour are sought is passed to the function via 
 parameter `int *order`.
 
 The return value is given via argument `int *vector`
 
 */


void c_balanced_clustering(
  double *data, int *N, int *M, 
  int *K, int *order, int *vector, int *mem_error) {
  
        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of variables per data point
        const size_t k = (size_t) *K; // number of clusters
        //const size_t n_per_c = n / k;
        // Some book-keeping variables to track memory error
        int mem_error_points = 0;
        
        // Set up array of data points, fill it, return if memory runs out
        struct cl_element POINTS[n];
        mem_error_points = set_up_list(data, n, m, POINTS);
        
        if (mem_error_points == 1) {
                *mem_error = 1;
                return;
        }
        
        // Set up linked list for data points
        struct double_node *HEAD = malloc(sizeof(struct double_node));
        if (HEAD == NULL) {
                *mem_error = 1;
                return;
        }
        HEAD->next = NULL;
        HEAD->prev = NULL;
        
        // 
        struct double_node *PTR_ARRAY[n];
        // Fill both data structures: (a) linked list and (b) array of pointers
        for (int i = (n-1); i >= 0; i--) { // traverse from the end
                struct double_node* new = insert_double_node(HEAD, &POINTS[order[i]]);
                if (new == NULL) {
                        *mem_error = 1;
                        //TODO: free memory
                        return;
                }
                // now add pointer to the new node to the pointer array. 
                // Careful indexing needed
                PTR_ARRAY[POINTS[order[i]].ID] = new;
        }
}


/* Insert a double node next to HEAD 
 * Returns the pointer to the new node so it can be stored in an array
 * 
 * Returns NULL if the memory allocation fails.
 */
struct double_node* insert_double_node(struct double_node *HEAD, struct cl_element *POINT) {
        struct double_node* new = malloc(sizeof(struct double_node));
        if (new == NULL) {
                return NULL;
        }
        new->prev = HEAD;
        new->next = HEAD->next;
        HEAD->next = new;
        new->data = POINT;
        //printf("%lf \n", new->data->values[0]);
        return new;
}


/* Extracted method that fill data points into array of struct `element`
 * param `double *data` pointer to original data array of length n
 * param `size_t m`: Number of data points
 * param `size_t m`: Number of variables per data point
 * param `struct element POINTS[n]`: Array to be filled with data points
 * 
 * return: `0` if all data points were successfully stored; `1` if not.
 * 
 */
int set_up_list(double *data, size_t n, size_t m, struct cl_element POINTS[n]) {
        // Create offset variable to correctly read out data points
        int m_ptr[m];
        for (size_t i = 0; i < m; i++) {
                m_ptr[i] = i * n;
        }
        
        // Size of a data vector per element:
        size_t data_size = m * sizeof(POINTS[0].values[0]);
        
        for (size_t i = 0; i < n; i++) {
                POINTS[i].cluster = 0; // init as 0, is written later
                POINTS[i].ID = i;
                POINTS[i].values = (double*) malloc(data_size);
                if (POINTS[i].values == NULL) {
                        //free_points(n, POINTS, i); // TODO: free for cl_element
                        return 1;
                } 
                // Fill data into `element`:
                for (size_t j = 0; j < m; j++) {
                        POINTS[i].values[j] = data[m_ptr[j]++];
                }
        }
        return 0;
}
 

/* Must implement the following functionality:
 * - remove item from linked list (is this even reasonably possible?! - I guess I need 
 *   another array as in the clustering functions to switch between list and array
 *   representation...). Well, no, this should only later be done when all neighbours 
 *   are definitely found and there are no more changes! Thus, I need a PTR_TO_NODES, 
 *   then it should easily work.
 * - add to the two arrays, in correct order! */
// insert_into_cluster(tmp, cluster_ids, cluster_distances);
