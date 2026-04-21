
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include "declarations.h"

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
        double *d = (double *) malloc(m * sizeof(double));     /* d = x - y */
        if (!d) { perror("malloc"); exit(EXIT_FAILURE); }
        double *t = (double *)malloc(m * sizeof(double));      /* t = Σ⁻¹ * d */
        if (!t) { perror("malloc"); exit(EXIT_FAILURE); }

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

        // clean up
        free(d);
        free(t);

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
