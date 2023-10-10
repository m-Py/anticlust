#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "declarations.h"

// Compute similarities a la clique partitioning 
void cpp_similarities_(int *data, int *N, int *M, int *output) {
        const size_t n = (size_t) *N;
        const size_t m = (size_t) *M;
        size_t counter = 0;
        for (size_t i = 0; i < n-1; i++) {
                for (size_t j = i+1; j < n; j++) {
                  int sum_disagreements = 0;
                  for (size_t z = 0; z < m; z++) {
                          sum_disagreements += data[one_dim_index(i, z, n)] != data[one_dim_index(j, z, n)];
                  }
                  output[counter] = m - 2 * sum_disagreements;
                  counter++;
                }
        }
}
