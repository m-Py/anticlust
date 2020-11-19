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
        const size_t n_per_c = n / k;
        // Some book-keeping variables to track memory error
        int mem_error_points = 0;
        int mem_error_cluster_lists = 0;
        int false = 0;
        
        // Set up array of data points, fill it, return if memory runs out
        struct element POINTS[n];
        mem_error_points = fill_data_points(
                data, n, m, POINTS, vector, // "vector" actually not needed, is hacky
                &false, vector
        );
        
        if (mem_error_points == 1) {
                *mem_error = 1;
                return;
        }
        
        // Set up linked list for data points
        struct node *HEAD = (struct node*) malloc(sizeof(struct node));
        if (HEAD == NULL) {
                *mem_error = 1;
                return;
        }
        HEAD->next = NULL;
        for (int i = (n-1); i >= 0; i--) { // traverse from the end
                struct node* tmp = append_to_cluster(HEAD, &POINTS[order[i]]);
                if (tmp == NULL) {
                        *mem_error = 1;
                        return;
                }
                printf("%lf \n", POINTS[order[i]].values[0]);
        }
        
        // Set up array of pointers-to-nodes, return if memory runs out
        struct node *PTR_NODES[n];
        mem_error_cluster_lists = fill_cluster_lists(
            n, k, clusters, 
            POINTS, PTR_NODES, CLUSTER_HEADS
        );
        
        // TODO: test if memory allocation failed and return if it did        
        
        // Iterate through list, find neighbours for each target element
        int cluster = 0; // counter for the clusters 
        while (HEAD->next != NULL) {
                // Use another list for the current cluster
                // 1. Pop next element from list
                struct node *target_item = HEAD->next; // for this, look for neighbours
                struct node *compare_to = target_item->next;
                double worst_distance = euclidean_squared(
                        target_item->data->values, compare_to->data->values, m
                );
                
                // Use two arrays to keep count of the neighbours
                // (a) one includes IDs
                // (b) one includes distances to the target item
                int cluster_ids[n_per_c]; 
                double cluster_distances[n_per_c];
                
                // Fill first two members to cluster, one may be replaced later
                cluster_ids[0] = target_item->data->ID;
                cluster_ids[1] = compare_to->data->ID;
                cluster_distances[0] = 0;
                cluster_distances[1] = worst_distance;
                target_item->data->cluster = cluster;
                compare_to->data->cluster = cluster; // may be reverted
                
                int members_in_cluster = 2; // counter for the members in a cluster
                
                struct node* tmp = compare_to->next;
                
                while (tmp->next != NULL) {
                        double tmp_distance = euclidean_squared(
                                target_item->data->values, tmp->data->values, m
                        );
                        // item is inserted if its closer to target item than the current
                        // worst item in the cluster, OR if the cluster is not yet filled
                        // entirely
                        if (tmp_distance < worst_distance || members_in_cluster < n_per_c) {
                                // to be implemented: insert into cluster logic
                                insert_into_cluster(HEAD, tmp, cluster_ids, cluster_distances);
                                if (tmp_distance > worst_distance) {
                                       worst_distance = tmp_distance;
                                }
                                if (members_in_cluster < n_per_c) {
                                        members_in_cluster++;
                                }
                        }
                        tmp = tmp->next;
                }
                
                // TODO: Add cluster affiliation to every data point
                
                /* TODO: Use `cluster_ids` array and `PTR_TO_NODES` array to remove the nodes 
                 * from the primary linked list */
                remove_from_linked_list() // to be implemented
                
                
                // at the end of this loop, HEAD->next must be updated
                // increment `cluster`
                cluster++;
                
        }
        
}


/* Must implement the following functionality:
 * - remove item from linked list (is this even reasonably possible?! - I guess I need 
 *   another array as in the clustering functions to switch between list and array
 *   representation...). Well, no, this should only later be done when all neighbours 
 *   are definitely found and there are no more changes! Thus, I need a PTR_TO_NODES, 
 *   then it should easily work.
 * - add to the two arrays, in correct order!
// insert_into_cluster(tmp, cluster_ids, cluster_distances);
