
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

// Function for k-means anticlustering using the Mahalanobis Distance
void fast_MD_anticlustering(double *data, int *N, int *M, int *K, int *frequencies,
                            int *clusters, int *partners, int *k_neighbours,
                            double *inv_cov) {

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
                OBJ_BY_CLUSTER[i] = mahalanobis_inv(OVERALL_CENTROID, CENTERS[i], m, inv_cov);
        }

        /* Some variables for bookkeeping during the optimization */
        size_t best_partner;
        size_t best_cluster;
        double tmp_center1[m];
        double tmp_center2[m];
        double tmp_obj_cl1;
        double tmp_obj_cl2;
        double best_obj_cl1;
        double best_obj_cl2;
        double best_centers_cl1[m];
        double best_centers_cl2[m];
        double tmp_reduction; // encode if the exchange led to reduction in objective (this is a minimization problem)

        /* Start main iteration loop for exchange procedure */
        size_t id_current_exch_partner = 0;

        /* 1. Level: Iterate through `n` data points */
        for (size_t i = 0; i < n; i++) {
                int cl1 = clusters[i];

                // Initialize `best` variable for the i'th item
                double best_reduction = 0;
                int exchange_cluster_found = 0;

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
                                copy_array(m, CENTERS[cl1], tmp_center1);
                                copy_array(m, CENTERS[cl2], tmp_center2);
                                tmp_obj_cl1 = OBJ_BY_CLUSTER[cl1];
                                tmp_obj_cl2 = OBJ_BY_CLUSTER[cl2];

                                fast_update_one_center(
                                        i, j, n, m,
                                        data, tmp_center1,
                                        frequencies[cl1]
                                );

                                fast_update_one_center(
                                        j, i, n, m,
                                        data, tmp_center2,
                                        frequencies[cl2]
                                );

                                // Update objective
                                tmp_obj_cl1 = mahalanobis_inv(OVERALL_CENTROID, tmp_center1, m, inv_cov); // should be smaller than before
                                tmp_obj_cl2 = mahalanobis_inv(OVERALL_CENTROID, tmp_center2, m, inv_cov);

                                // Update objective
                                tmp_reduction = tmp_obj_cl1 * frequencies[cl1] + tmp_obj_cl2 * frequencies[cl2] -
                                OBJ_BY_CLUSTER[cl1] * frequencies[cl1] - OBJ_BY_CLUSTER[cl2] * frequencies[cl2];

                                // Update `best` variables if objective was improved
                                if (tmp_reduction < best_reduction) {
                                        best_obj_cl1 = tmp_obj_cl1;
                                        best_obj_cl2 = tmp_obj_cl2;
                                        copy_array(m, tmp_center1, best_centers_cl1);
                                        copy_array(m, tmp_center2, best_centers_cl2);
                                        best_partner = j;
                                        best_cluster = cl2;
                                        exchange_cluster_found = 1;
                                        best_reduction = tmp_reduction;
                                }
                        }
                }

                // Only if objective is improved: Do the swap
                if (exchange_cluster_found) {
                        fast_swap(clusters, i, best_partner);
                        // Update the "global" variables
                        OBJ_BY_CLUSTER[cl1] = best_obj_cl1;
                        OBJ_BY_CLUSTER[best_cluster] = best_obj_cl2;
                        copy_array(m, best_centers_cl1, CENTERS[cl1]);
                        copy_array(m, best_centers_cl2, CENTERS[best_cluster]);
                        // here for local-maximum method: set flag that an improvement occurred!
                }

        }
        return;
}

/* Squared Mahalanobis Distance between two arrays of same length
 *
 * param *x: Array / pointer to first element
 * param *y: Array / pointer to second element
 * param m: length of the two arrays
 * param *inv_cov: inverted covariance matrix (is computed in R), in column‑major order
 *
 * return: The squared MD distance distance
 *
 */
double mahalanobis_inv(double *x,
                       double *y,
                       size_t m,
                       double *inv_cov) {


        // Vectors to track computation of MD
        double d[m];     /* d = x - y */
        double t[m];      /* t = Σ⁻¹ * d */

        for (size_t i = 0; i < m; ++i) {
                d[i] = x[i] - y[i];
        }

        // First step of MD: matrix vector produkt
        for (size_t i = 0; i < m; ++i) {
                double sum = 0.0;
                for (size_t j = 0; j < m; ++j)
                        sum += inv_cov[i + j * m] * d[j];   /* Σ⁻¹_ij * d_j */
                t[i] = sum;
        }

        // Second step of MD: scalar product
        double dot = 0.0;
        for (size_t i = 0; i < m; ++i) {
                dot += d[i] * t[i];
        }

        return dot;
}

// for directly calling from R
void mahalanobis_out(double *x,
                     double *y,
                     int *m,
                     double *inv_cov,
                     double *MD) {
        size_t M = (size_t) *m;
        *MD = mahalanobis_inv(x, y, M, inv_cov);
}
