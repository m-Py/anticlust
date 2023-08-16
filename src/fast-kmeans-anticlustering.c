
#include <math.h>
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
 * param *mem_error: This is passed with value 0 and only receives the value 1 
 *       if a memory error occurs when executing this function. The caller needs
 *       to test if this value is 1 after execution.
 * 
 * The return value is assigned to the argument `clusters`, via pointer
 * 
 * 
 * ============================ Some explanations ============================
 * 
 * This is a slightly different version of kmeans_anticlustering(), which does not
 * use categorical restrictions but has exchange partners for each element (which must be
 * passed to the function and are not generated here)
 * ===========================================================================
*/

void fast_kmeans_anticlustering(double *data, int *N, int *M, int *K, int *frequencies,
        int *clusters, int *partners, int *k_neighbours, int *mem_error) {
        
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
        const size_t kn = (size_t) *k_neighbours; // number of clusters
        
        // Do the same for the original data points
        double data_pts[n][m];
        int m_ptr[m];
        
        // Column offsets (to convert one-dimensional array to Row/Col major)
        for(int i = 0; i < m; i++) {
                m_ptr[i] = i * n;
        }
        
        // Reconstruct the data points as N x M matrix
        for(int i = 0; i < n; i++) {
                for(int j = 0; j < m; j++) {
                        data_pts[i][j] = data[m_ptr[j]++];
                }
        }

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
        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < m; j++) {
            CENTERS[clusters[i]][j] += data_pts[i][j];
          }
        }
        // To get cluster centers: Divide by number of elements 
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
        for (size_t i = 0; i < n; i++) {
          OBJ_BY_CLUSTER[clusters[i]] += euclidean_squared(data_pts[i], CENTERS[clusters[i]], m);
        }

        /* SUM OF SQUARES ACROSS ALL CLUSTERS */
        double SUM_OBJECTIVE = 0;
        for (size_t i = 0; i < k; i++) {
          SUM_OBJECTIVE += OBJ_BY_CLUSTER[i];
        }
        printf("Variance is %g \n", SUM_OBJECTIVE);

        /* Some variables for bookkeeping during the optimization */
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
          double best_obj = 0;
          copy_matrix(k, m, CENTERS, best_centers);
          copy_array(k, OBJ_BY_CLUSTER, best_objs);
          size_t j;
          /* 2. Level: Iterate through the exchange partners */
          for (size_t u = 0; u < n; u++) {
            // Get index of current exchange partner
            j = u;
            int cl2 = clusters[j];
            // no swapping attempt if in the same cluster:
            if (cl1 == cl2) { 
              continue;
            }
            
            // Initialize `tmp` variables for the exchange partner:
            copy_matrix(k, m, CENTERS, tmp_centers);

            fast_update_centers(
              i, j, n, m, k,
              data_pts, cl1, cl2,
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
              // do not compute if clusters[f] is cl1 or cl2; optimization for later
              tmp_objs[clusters[f]] += euclidean_squared(data_pts[f], tmp_centers[clusters[f]], m);
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
          id_current_exch_partner++;
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