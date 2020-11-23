#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
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
        const size_t n_per_c = n / k;
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
        struct double_node *HEAD = malloc(sizeof(struct double_node*));
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
        

        
        // Iterate through list, find neighbours for each target element
        int cluster = 0; // counter for the clusters
        struct cl_node *CL_HEAD = malloc(sizeof(struct cl_node*));
        while (HEAD->next != NULL) {
                // Use another list for the current cluster
                // 1. Pop next element from list
                struct double_node *TARGET = HEAD->next; // for this, look for neighbours
                double worst_distance = 0; 
                
                // Set up linked list representing this cluster
                // TODO: Check memory
                CL_HEAD->element = TARGET;
                CL_HEAD->next = NULL;
                CL_HEAD->distance = 0;
                
                // Fill first two members to cluster, one may be replaced later
                int members_in_cluster = 1; // counter for the members in a cluster
                
                // Iterate through the double linked list, containing all (remaining)
                // elements. Append to cluster list if new element is a nearest neighbour
                struct double_node* tmp = TARGET->next;
                while (tmp != NULL) {
                        double tmp_distance = euclidean_squared(
                                TARGET->data->values, tmp->data->values, m
                        );
                        // item is inserted if its closer to target item than the current
                        // worst item in the cluster, OR if the cluster is not yet filled
                        // entirely
                        if (tmp_distance < worst_distance || members_in_cluster < n_per_c) {
                                insert_into_cluster(
                                        CL_HEAD,
                                        tmp, 
                                        tmp_distance,
                                        cluster
                                );
                                if (tmp_distance > worst_distance) {
                                       worst_distance = tmp_distance;
                                }
                                members_in_cluster++;
                        }
                        tmp = tmp->next;
                }
                
                // Iterate through cluster list and remove the assigned elements
                // from the item pool, so they are no longer used when searching 
                // for nearest neighbours. Free the cluster list, a new one is 
                // in the next iteration of the outer loop.
                struct cl_node *tmp2 = CL_HEAD;
                size_t counter = 0;
                while (tmp2 != NULL) {
                        size_t index = tmp2->element->data->ID;
                        struct double_node *current = PTR_ARRAY[index];
                        if (counter < n_per_c) {
                                // Update cluster of the current element
                                current->data->cluster = cluster;
                                // delete from doubled linked list
                                // Case 1: current->next is NULL (i.e., end of list)
                                if (current->next == NULL) {
                                        current->prev->next = NULL;
                                } else {
                                        current->prev->next = current->next;
                                        current->next->prev = current->prev;
                                }
                        }
                        // redefine distance as INFINITY so node can be reused
                        tmp2->distance = INFINITY; 
                        //struct cl_node *del2 = tmp2;
                        tmp2 = tmp2->next;
                        //free(del2); // can be freed
                        counter++;
                }
                cluster++;
        } 
        // Write output
        for (size_t i = 0; i < n; i++) {
                vector[i] = PTR_ARRAY[i]->data->cluster;
        }
        
}

/* Must implement the following functionality:
 * - add to the two arrays, in correct order! */
void insert_into_cluster(
        struct cl_node *HEAD, 
        struct double_node *node,
        double distance,
        int cluster
) {
        struct cl_node *tmp = HEAD;
        while (tmp != NULL) {
                // Just insert if we are at the end of the list
                if (tmp->next == NULL) {
                      // Make new cl_node for new element
                      struct cl_node *new = malloc(sizeof(struct cl_node*));
                      //TODO: memory check
                      new->distance = distance;
                      new->element = node;
                      tmp->next = new;
                      tmp->next->next = NULL;
                      return; 
                }
                // reuse node if possible
                if (tmp->next->distance == INFINITY) {
                        tmp->next->distance = distance;
                        tmp->next->element = node;
                        return;
                }
                // here: insert into the list, somewhere in the middle
                // if new distance is lower than next distance, insert here
                if (distance < tmp->next->distance) {
                        struct cl_node *new = malloc(sizeof(struct cl_node*));
                        new->distance = distance;
                        new->element = node;
                        struct cl_node *safe = tmp->next;
                        tmp->next = new;
                        new->next = safe;
                        return;
                }
                tmp = tmp->next;
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
        if (new->next != NULL) {
                new->next->prev = new;
        }
        HEAD->next = new;
        new->data = POINT;
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
 
