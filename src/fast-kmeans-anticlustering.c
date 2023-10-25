
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

/* Exchange Method for K-Means Anticlustering
 * param *data: vector of data points (in R, this is a data frame / matrix,
 *         but it is passed as a pointer array (vector) to C; here, we use
 *         "clever" indexing on the one-dimensional vector to address the correct 
 *         data)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *M: The number of features (i.e., number of "cols" in *data)
 * param *K: The number of clusters
 * param *frequencies: The number of elements per cluster, i.e., an array
 *         of length *K.
 * param *clusters: An initial assignment of elements to clusters,
 *         array of length *N (has to consists of integers between 0 and (K-1) 
 *         - this has to be guaranteed by the caller)
 * param *partners: A pointer array of length N * k_neighbours, indicating for each 
 *        element which other elements are exchange partners. The first `k_neighbours`
 *        entries are the exchange partners of element 1, the next belong to element
 *        2, and so forth. If an entry is N, the element is skipped (only indices up to N-1 work, 
 *        so N is used to indicate that no exchange partners follow).
 * param *k_neighbours: The number of exchange partners per element.
 * 
 * The return value is assigned to the argument `clusters`, via pointer
 * 
 * 
 * ============================ Some explanations ============================
 * 
 * This is a new implementation of k-means anticlustering that uses an alternative
 * computation of the k-means objective: Sum of squared Euclidean 
 * distances between cluster centers and overall centroid; it should
 * be minimized.
 * ===========================================================================
*/

void fast_kmeans_anticlustering(double *data, int *N, int *M, int *K, int *frequencies,
        int *clusters, int *partners, int *k_neighbours) {
        
        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of variables per data point
        const size_t k = (size_t) *K; // number of clusters
        const size_t kn = (size_t) *k_neighbours; // number of clusters
        
        /* INITIALIZE OBJECTIVE */
        // 1 OVERALL CENTROID
        double OVERALL_CENTROID[m];
        init_overall_centroid(m, n, OVERALL_CENTROID, data);
        
        /* K CLUSTER CENTERS */
        // Set up matrix of cluster centers
        double CENTERS[k][m];
        init_centers(k, m, n, CENTERS, clusters, frequencies, data);
        
        /* DISTANCES BETWEEN CLUSTER CENTERS AND OVERALL CENTROID */
        double OBJ_BY_CLUSTER[k]; 
        for (size_t i = 0; i < k; i++) {
          OBJ_BY_CLUSTER[i] = euclidean_squared(OVERALL_CENTROID, CENTERS[i], m);
        }
        /* K-means objective in this setting: Sum of squared Euclidean 
         * distances between cluster centers and overall centroid */
        double SUM_OBJECTIVE = weighted_array_sum(k, frequencies, OBJ_BY_CLUSTER); 

        /* Some variables for bookkeeping during the optimization */
        size_t best_partner;
        double tmp_centers[k][m];
        double best_centers[k][m];
        double tmp_objs[k];
        double best_objs[k];
        double tmp_obj;
        
        /* Start main iteration loop for exchange procedure */
        size_t id_current_exch_partner = 0;
        
        
        /* 1. Level: Iterate through `n` data points */
        for (size_t i = 0; i < n; i++) {
          int cl1 = clusters[i];
          
          // Initialize `best` variable for the i'th item
          double best_obj = SUM_OBJECTIVE + 1; // we have a minimization problem; init as large value
          copy_matrix(k, m, CENTERS, best_centers);
          copy_array(k, OBJ_BY_CLUSTER, best_objs);
          
          size_t j;
          /* 2. Level: Iterate through the exchange partners */
          
          for (size_t u = 0; u < kn; u++) {
            // Get index of current exchange partner
            j = partners[id_current_exch_partner];
            id_current_exch_partner++; // this just counts upwards across all exchange partners, ignores matrix-like structure
            if (j != n) { // no exchange partners any more
              int cl2 = clusters[j];
              // no swapping attempt if in the same cluster:
              if (cl1 == cl2) {
                continue;
              }
              
              // Initialize `tmp` variables for the exchange partner:
              // Center matrix, is updated after swap
              // Initialize `tmp` variables for the exchange partner:
              copy_matrix(k, m, CENTERS, tmp_centers);
              copy_array(k, OBJ_BY_CLUSTER, tmp_objs);
              
              fast_update_centers(
                i, j, n, m, k,
                data, cl1, cl2,
                tmp_centers, 
                frequencies
              );
              
              // Update objective
              tmp_objs[cl1] = euclidean_squared(OVERALL_CENTROID, tmp_centers[cl1], m);
              tmp_objs[cl2] = euclidean_squared(OVERALL_CENTROID, tmp_centers[cl2], m);
              tmp_obj = weighted_array_sum(k, frequencies, tmp_objs);
              
              // Update `best` variables if objective was improved
              if (tmp_obj < best_obj) {
                best_obj = tmp_obj;
                copy_matrix(k, m, tmp_centers, best_centers);
                copy_array(k, tmp_objs, best_objs);
                best_partner = j;
              }
            }
          }
          
          // Only if objective is improved: Do the swap
          if (best_obj < SUM_OBJECTIVE) {
            fast_swap(clusters, i, best_partner);
            // Update the "global" variables
            SUM_OBJECTIVE = best_obj;
            copy_matrix(k, m, best_centers, CENTERS);
            copy_array(k, best_objs, OBJ_BY_CLUSTER);
            // here for local-maximum method: set flag that an improvement occurred!
          }
          
        }
        return;
}


/* Update cluster centers for simpler implementation not using cluster lists */
void fast_update_centers(size_t i, size_t j, size_t n, size_t m, size_t k, double *data, 
                         int cl1, int cl2, double CENTERS[k][m], int *frequencies) {
  
  for (int u = 0; u < m; u++) {
    double added_to_cl1 = data[one_dim_index(j, u, n)] / frequencies[cl1];
    double removed_from_cl1 = data[one_dim_index(i, u, n)] / frequencies[cl1];
    
    double added_to_cl2 = data[one_dim_index(i, u, n)] / frequencies[cl2];
    double removed_from_cl2 = data[one_dim_index(j, u, n)] / frequencies[cl2];
    
    // Update first cluster center
    CENTERS[cl1][u] = CENTERS[cl1][u] + added_to_cl1;
    CENTERS[cl1][u] = CENTERS[cl1][u] - removed_from_cl1;
    // Update second cluster center
    CENTERS[cl2][u] = CENTERS[cl2][u] - removed_from_cl2;
    CENTERS[cl2][u] = CENTERS[cl2][u] + added_to_cl2;
  }
}


void init_overall_centroid(size_t m, size_t n, double OVERALL_CENTROID[m], double* data) {
  // Init as Zero
  for (size_t j = 0; j < m; j++) {
    OVERALL_CENTROID[j] = 0; 
  }
  /* NOW JUST SUM UP ALL VALUES */
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      OVERALL_CENTROID[j] += data[one_dim_index(i, j, n)];  
    }
  }
  /* DIVIDE BY N */
  for (size_t j = 0; j < m; j++) {
    OVERALL_CENTROID[j] = OVERALL_CENTROID[j] / n;  
  }
}

  
  
void init_centers(size_t k, size_t m, size_t n, double CENTERS[k][m], int* clusters, int* frequencies, double* data) {
  // Init as Zero
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < m; j++) {
      CENTERS[i][j] = 0; 
    }
  }
  /* SUM UP FEATURE VALUES BY CLUSTER */
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      // data[n * j + i] == data[i, j] in the original matrix, but in C we only have a 1-dim pointer
      CENTERS[clusters[i]][j] += data[one_dim_index(i, j, n)];  
    }
  }
  /* To get cluster centers: Divide by number of elements */
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < m; j++) {
      CENTERS[i][j] = CENTERS[i][j] / frequencies[i];
    }
  }
}

void fast_swap(int *clusters, size_t i, size_t j) {
  int tmp = clusters[i];
  clusters[i] = clusters[j];
  clusters[j] = tmp;
}

size_t one_dim_index(size_t i, size_t j, size_t n) {
  return n * j + i;
}

