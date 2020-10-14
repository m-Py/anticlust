
#include "declarations.h"
#include <stdlib.h> 
#include <math.h>

/* Compute sum of squared errors between cluster center and data points (i.e.
 * variance) for each of the k clusters. The input array `VAR_OBJECTIVE` is 
 * changed through this function. 
 */
void objective_by_cluster(size_t m, size_t k, double OBJ_BY_CLUSTER[k], 
                          double CENTERS[k][m], struct node *HEADS[k]) {
        for (size_t i = 0; i < k; i++) {
                OBJ_BY_CLUSTER[i] = cluster_var(m, HEADS[i], CENTERS[i]);
        }
}

/* Compute variance for a cluster
 * param `size_t m`: Number of variables per data point
 * param `struct node *HEAD`: Pointer to a cluster list HEAD
 * param `double center[m]`: Array of mean feature values in the cluster
 */
double cluster_var(size_t m, struct node *HEAD, double center[m]) {
        double sum = 0;
        struct node *tmp = HEAD->next;
        while (tmp != NULL) {
                sum += euclidean_squared(center, tmp->data->values, m);
                tmp = tmp->next;
        }
        return sum;
}


/* Compute cluster center for one cluster
 * 
 * param `size_t m`: Number of variables per data point
 * param `double center[m]`: Empty array of cluster centers
 * param `struct node *HEAD`: Pointer to a cluster list HEAD
 * param `int freq`: Number of elements in the cluster
 * 
 * The input array `center` is changed through this function to represent
 * one cluster center.
 * 
 */
void compute_center(size_t m, double center[m], struct node *HEAD, int frequency) {
        // Initialize center matrix as 0:
        for (size_t i = 0; i < m; i++) {
                center[i] = 0; 
        }
        
        struct node *tmp = HEAD->next; 
        while (tmp != NULL) {
                for (size_t i = 0; i < m; i++) {
                        center[i] = center[i] + tmp->data->values[i];
                }
                tmp = tmp->next;
        } 
        // To get cluster centers: Divide by number of elements 
        for (size_t i = 0; i < m; i++) {
                center[i] = center[i] / frequency;
        }
}


/* Squared Euclidean Distance between two arrays of same length
 * 
 * param *x: Array / pointer to first element
 * param *y: Array / pointer to second element
 * param m: length of the two arrays
 * 
 * return: The squared euclidean distance
 * 
 */

double euclidean_squared(double *x, double *y, size_t m) {
        double sum = 0;
        for (size_t i = 0; i < m; i++) {
                sum = sum + pow(x[i] - y[i], 2);
        }
        return sum;
}


/* Update cluster centers after a swap between two nodes in two cluster lists */
void update_centers(size_t k, size_t m, double CENTERS[k][m],
                    struct node *one, struct node *two, int *frequencies) {
        size_t cl1 = one->data->cluster;
        size_t cl2 = two->data->cluster;
        for (int i = 0; i < m; i++) {
                double added_to_cl1 = two->data->values[i] / frequencies[cl1];
                double removed_from_cl1 = one->data->values[i] / frequencies[cl1];
                
                double added_to_cl2 = one->data->values[i] / frequencies[cl2];
                double removed_from_cl2 = two->data->values[i] / frequencies[cl2];
                
                // Update first cluster center
                CENTERS[cl1][i] = CENTERS[cl1][i] + added_to_cl1;
                CENTERS[cl1][i] = CENTERS[cl1][i] - removed_from_cl1;
                // Update second cluster center
                CENTERS[cl2][i] = CENTERS[cl2][i] - removed_from_cl2;
                CENTERS[cl2][i] = CENTERS[cl2][i] + added_to_cl2;
        }
}
