#include <stdlib.h> 


/* Copy one array into another */
void copy_array(size_t n, double origin[n], double target[n]) {
        for (int i = 0; i < n; i++) {
                target[i] = origin[i];
        }
}

/* Copy one matrix into another */
void copy_matrix(size_t n, size_t m, double origin[n][m], double target[n][m]) {
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                        target[i][j] = origin[i][j];
                }
        }
}

/* Compute the sum of an array */
double array_sum(size_t k, double ARRAY[k]) {
        double sum = 0;
        for (size_t i = 0; i < k; i++) {
                sum += ARRAY[i];
        }
        return sum;
}
 
