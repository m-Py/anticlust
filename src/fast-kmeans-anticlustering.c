
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

/* Exchange Method for K-Means Anticlustering
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
 * This is a slightly different version of kmeans_anticlustering(), which does not
 * use categorical restrictions but has exchange partners for each element (which must be
 * passed to the function and are not generated here).
 * 
 * The code is horrible because it tries to reduce additional function calls, 
 * which seemed to terminate my R sessions for really large data (and it should be
 * usable for really large data sets). Probably it was not the function calls,
 * but me writing additional (and not necessarily needed data), but here we are.
 * Now it works, even for N > 250000 (where the old C implementation failed), 
 * and I am quite unwilling to change a lot.
 * ===========================================================================
*/

void fast_kmeans_anticlustering(double *data, int *N, int *M, int *K, int *frequencies,
        int *clusters, int *partners, int *k_neighbours) {
        
        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of variables per data point
        const size_t k = (size_t) *K; // number of clusters
        const size_t kn = (size_t) *k_neighbours; // number of clusters
        
        /* INITIALIZE OBJECTIVE */
        
        /* CLUSTER CENTERS */
        // Set up matrix of cluster centers
        double CENTERS[k][m];
        // Init as Zero // It seems I cannot use functions here due to Stack Overflow for large N, so all is in the main function...
        for (size_t i = 0; i < k; i++) {
          for (size_t j = 0; j < m; j++) {
            CENTERS[i][j] = 0; 
          }
        }
        /* SUM UP FEATURE VALUES BY CLUSTER */
        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < m; j++) {
            // data[n * j + i] == data[i, j] in the original matrix, but in C we only have a 1-dim pointer
            CENTERS[clusters[i]][j] += data[n * j + i];  
          }
        }
        /* To get cluster centers: Divide by number of elements */
        for (size_t i = 0; i < k; i++) { // excuse my reuse of indices, I hate this
          for (size_t j = 0; j < m; j++) {
            CENTERS[i][j] = CENTERS[i][j] / frequencies[i];
          }
        }

        /* SUM OF SQUARES BY CLUSTER */
        double OBJ_BY_CLUSTER[k]; 
        for (size_t i = 0; i < k; i++) {
          OBJ_BY_CLUSTER[i] = 0; // init as zero
        }
        /* Compute Squared Euclidean Distance between each data point and its center -> sum up */
        for (size_t i = 0; i < n; i++) {
          double distance_squared = 0;
          for (size_t j = 0; j < m; j++) {
            double diff = data[n * j + i] - CENTERS[clusters[i]][j];
            distance_squared += diff * diff;
          }
          OBJ_BY_CLUSTER[clusters[i]] += distance_squared;
        }
        

        /* SUM OF SQUARES ACROSS ALL CLUSTERS */
        double SUM_OBJECTIVE = 0;
        for (size_t i = 0; i < k; i++) {
          SUM_OBJECTIVE += OBJ_BY_CLUSTER[i];
        }
        

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
          
          // Initialize `best` variables (alas not by function call for really large data)
          double best_obj = 0;
          // Best Center matrix
          for (int i2 = 0; i2 < k; i2++) {
            for (int j2 = 0; j2 < m; j2++) {
              best_centers[i2][j2] = CENTERS[i2][j2];
            }
          }
          // Best objectives by cluster
          for (int i3 = 0; i3 < k; i3++) { // WTF indices...
            best_objs[i3] = OBJ_BY_CLUSTER[i3];
          }
          
          size_t j;
          /* 2. Level: Iterate through the exchange partners */

          for (size_t u = 0; u < kn; u++) {
            // Get index of current exchange partner
            j = partners[id_current_exch_partner];
            if (j == n) { // no exchange partners any more
              continue;
            }
            id_current_exch_partner++; // this just counts upwards across all exchange partners, ignores matrix-like structure
            
            int cl2 = clusters[j];
            // no swapping attempt if in the same cluster:
            if (cl1 == cl2) {
              continue;
            }
            
            // Initialize `tmp` variables for the exchange partner:
            // Center matrix, is updated after swap
            for (int i3 = 0; i3 < k; i3++) {
              for (int j3 = 0; j3 < m; j3++) {
                tmp_centers[i3][j3] = CENTERS[i3][j3];
              }
            }

            fast_update_centers(
              i, j, n, m, k,
              data, cl1, cl2,
              tmp_centers, 
              frequencies
            );

            fast_swap(clusters, i, j);
            // Update objective
            /* SUM OF SQUARES BY CLUSTER */
            for (size_t f = 0; f < k; f++) {
              tmp_objs[f] = 0; // init as zero
            }
            for (size_t f = 0; f < n; f++) {
              double distance_squared = 0;
              for (size_t f2 = 0; f2 < m; f2++) {
                double diff = data[one_dim_index(f, f2, n)] - tmp_centers[clusters[f]][f2];
                distance_squared += diff * diff;
              }
              tmp_objs[clusters[f]] += distance_squared;
            }

            tmp_obj = array_sum(k, tmp_objs);
            
            // Update `best` variables if objective was improved
            if (tmp_obj > best_obj) {
              best_obj = tmp_obj;
              copy_matrix(k, m, tmp_centers, best_centers);
              copy_array(k, tmp_objs, best_objs);
              best_partner = j;
            }
            
            // Swap back to test next exchange partner
            fast_swap(clusters, i, j);
          }
          
          // Only if objective is improved: Do the swap
          if (best_obj > SUM_OBJECTIVE) {
            fast_swap(clusters, i, best_partner);
            // Update the "global" variables
            SUM_OBJECTIVE = best_obj;
            copy_matrix(k, m, best_centers, CENTERS);
            copy_array(k, best_objs, OBJ_BY_CLUSTER);
          }
        }
        return;
}


void fast_swap(int *clusters, size_t i, size_t j) {
  int tmp = clusters[i];
  clusters[i] = clusters[j];
  clusters[j] = tmp;
}

/* Function to print a matrix
 * @param N Number of rows in `matrix`
 * @param M Number of columns in `matrix`
 * @param matrix The matrix to be transposed, has dimension N x M
 * @return VOID
 */ 
void print_matrix(size_t N, size_t M, double matrix[N][M]) {
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < M; j++) {
      printf(" %10g ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("number of rows: %zu \nnumber of columns: %zu \n\n",
         N, M);
}

size_t one_dim_index(size_t i, size_t j, size_t n) {
  return n * j + i;
}

