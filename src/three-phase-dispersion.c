#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include "three-phase-header.h"

extern int beta_max;
extern int N, K;  // node number and group number
extern double** Distances;   // distance matrix between elements
extern double** DistancesT;
extern int* LB; // Lower bound of the number of elements in the group i 
extern int* UB; // Upper bound of the number of elements in the group i
extern double theta, theta_max, theta_min; 
extern double alpha;
extern int beta_min; 
extern int eta_max;

// main processure
extern int maxNumberIterations;
extern double start_time, Time_limit;
extern Solution S_best; //best solution
Solution *S_D; //S_i
Solution *O_D; //O_i in crossover

//double neighboorhood local search
extern double f_objective;
double *min_distance_per_cluster;
int **min_distance_tuple;
int* tuple1;
int* tuple2;

// Crosssover
double *min_distance_per_cluster_s1;
int **min_distance_tuple_s1;
double *min_distance_per_cluster_s2;
int **min_distance_tuple_s2;
extern int *vectorElement;
extern int *groupElement;
extern int *SelectGroup;
extern int *SelectEle;
extern int *tmpEle;
extern int *s1;
extern int *s2;
extern int *LBGroup;
extern int *UBGroup;
extern int *BigThanLB;
extern int *tmpUB;

//directed pertubation
extern int* Rd, * UnderLB; //Rd=R

void adding(int new_ind, int g, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void removing(int removed_ind, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void swapping(int ind1, int ind2, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
double evaluate_objective(double *s_min_distance_per_cluster);
void fill_arrays(int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void initialize_arrays(int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void recalculate_cluster_distance(int k, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void DoubleNeighborhoodLocalSearchDispersion(int s[], int SizeGroup[], double* objective);
void SearchAlgorithmDisperion(void);
void CrossoverDispersion(int partition1[], int partition2[], int score[], int scSizeGroup[]);
void DirectPerturbationDispersion(int eta_max, int s[], int SizeG[]);
void AssignMemoryDispersion(void);
void ReleaseMemoryDispersion(void);

/* TPSPD for Anticlustering Based on a Distance matrix
 * 
 * param *distannces: vector of data points (in R, this is a distance matrix,
 *         the matrix structure must be restored in C)
 * param *N_in: The number of elements (i.e., number of "rows" in *data)
 * param *K_in: The number of clusters. When lower_bound and upper_boound are set to the number of K, 
 *              the clusters will be equally sized. 
 * param *number_of_iterations: A number that defines how many times the steps in the search algorithm are repeated.
 * param *clusters: A predefined vector of length K specifies the number of elements in each cluster.
 *               If a default vector [-1] is provided, cluster sizes will be determined based 
 *               on the lower and upper bounds. When a cluster size array is provided, 
 *               the lower and upper bounds are ignored as they become redundant.
 * param *lower_bound: Minimum number of elements in each anticluster. 
 * param *upper_bound: Maximum number of elements in each anticluster.
 * param *Beta_max: The algorithm begins with a pool of random initial solutions of size beta_max. 
 *                   Over time, the size of the solution pool decreases linearly until it reaches beta_min.
 * param *elapsed_time: Measures the runtime of the algotihm (in seconds)
 * param *Theta_max: Parameter for the strength of undirected perturbation,
 *                   which decreases linearly over time from theta_max to theta_min..
 * param *Theta_min: Parameter for the strength of undirected perturbation, 
 *                   which decreases linearly over time from theta_max to theta_min..
 * param *Beta_min: The minimum solution pool size the algorithm should reach before making a determination.
 * param *Eta_max: A parameter that specifies how many times the steps in the direct perturbation are executed.
 * param *Alpha: Parameter for weitghing the discrimitation of a slighlty worse local optiomal child solution
 *               in Yang et al. set to 0.05 (might differ due to different implemetnation of calculation).
 * param *result: Calculated assignment of elements to clusters. Emptz vector.
 * param *objective: Value of objective function.
 * param *mem_error: This is passed with value 0 and only receives the value 1 
 *       if a memory error occurs when executing this function. The caller needs
 *       to test if this value is 1 after execution.
 * 
 * 
 * The return value is assigned to the argument `result`, via pointer
*/
void three_phase_search_dispersion(
                      double *distances, 
                      int *N_in,
                      int *K_in, 
                      int *number_of_iterations,
                      int *clusters,
                      int *upper_bound, 
                      int *lower_bound,
											int *Beta_max, 
											int *elapsed_time,
											double *Theta_max,
											double *Theta_min,
											int *Beta_min,
											int *Eta_max,
                                            double *Alpha,
											int *result,
											double *score,
											int *mem_error) {

  N = *N_in;
  K = *K_in;
  beta_max = *Beta_max;  
  theta = *Theta_max;
  theta_max = *Theta_max;
  beta_min = *Beta_min;
  eta_max = *Eta_max;
  Time_limit =  *elapsed_time;
  alpha = *Alpha;
  maxNumberIterations = *number_of_iterations;
  
  // Allocate memory for Distances and DistancesT arrays
  Distances = (double**)malloc(N * sizeof(double*));
  if (Distances == NULL) { *mem_error = 1; return; }
  DistancesT = (double**)malloc(N * sizeof(double*));
  if (DistancesT == NULL) { *mem_error = 1; return; }
  for (int i = 0; i < N; i++) {
    Distances[i] = (double*)malloc(N * sizeof(double));
    if (Distances[i] == NULL) { *mem_error = 1; return; }
    DistancesT[i] = (double*)malloc(N * sizeof(double));
    if (DistancesT[i] == NULL) { *mem_error = 1; return; }
  }

  //Allocates memory for DoubleSearhNeigboorhood
    tuple1 = (int *)malloc(2 * sizeof(int));
    tuple2 = (int *)malloc(2 * sizeof(int));
    
  // Fill Distances and DistancesT with values from input
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Distances[i][j] = distances[i * N + j];
      DistancesT[i][j] = 2 * distances[i * N + j];
    }
  }
  
  min_distance_per_cluster = (double *)malloc(K * sizeof(double));
  min_distance_tuple = malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    min_distance_tuple[i] = malloc(2 * sizeof(int));
  }
  min_distance_per_cluster_s1 = (double *)malloc(K * sizeof(double));
  min_distance_tuple_s1 = malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    min_distance_tuple_s1[i] = malloc(2 * sizeof(int)); // IMPORTANT: typo: missing s1 / s2 below
  }
  min_distance_per_cluster_s2 = (double *)malloc(K * sizeof(double));
  min_distance_tuple_s2 = malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    min_distance_tuple_s2[i] = malloc(2 * sizeof(int));
    
  }
  
  if (clusters[0] == -1 ) {
    // Allocate memory for LB and UB arrays
    LB = (int*)malloc(K * sizeof(int));
    if (LB == NULL) { *mem_error = 1; return; }
    UB = (int*)malloc(K * sizeof(int));
    if (UB == NULL) { *mem_error = 1; return; }
    for (int i = 0; i < K; i++) {
        LB[i] = *lower_bound;  // Assuming lower_bound is a pointer to an int
        UB[i] = *upper_bound;  // Assuming upper_bound is a pointer to an int
    }
  } else {
    LB = (int*)malloc(K * sizeof(int));
    if (LB == NULL) { *mem_error = 1; return; }
    UB = (int*)malloc(K * sizeof(int));
    if (UB == NULL) { *mem_error = 1; return; }
    for (int i = 0; i < K; i++) {
        LB[i] = clusters[i];  // Assuming lower_bound is a pointer to an int
        UB[i] = clusters[i];  // Assuming upper_bound is a pointer to an int
    }
  }

  AssignMemoryDispersion();
  if (*mem_error == 1) {
    return;
  }
  
  SearchAlgorithmDisperion();
  
  //save S_best -> solution with result
  for (int i = 0; i < N; i++){
    result[i] = S_best.s[i];
  }
  *score = S_best.objective;
  
  // Remember to free the allocated memory after use
  for (int i = 0; i < K; i++){
        free(min_distance_tuple[i]); min_distance_tuple[i] = NULL;
        free(min_distance_tuple_s1[i]); min_distance_tuple_s1[i] = NULL;
        free(min_distance_tuple_s2[i]); min_distance_tuple_s2[i] = NULL;
  }
  free(min_distance_per_cluster); min_distance_per_cluster = NULL;
  free(min_distance_tuple); min_distance_tuple = NULL;
  free(min_distance_per_cluster_s1); min_distance_per_cluster_s1 = NULL;
  free(min_distance_tuple_s1); min_distance_tuple_s1 = NULL;
  free(min_distance_per_cluster_s2); min_distance_per_cluster_s2 = NULL;
  free(min_distance_tuple_s2); min_distance_tuple_s2 = NULL;
  for (int i = 0; i < N; i++) {
    free(Distances[i]); Distances[i] = NULL;
    free(DistancesT[i]); DistancesT[i] = NULL;
  }
  free(Distances); Distances = NULL;
  free(DistancesT); DistancesT = NULL;
  free(LB); LB = NULL;
  free(UB); UB = NULL;
  free(tuple1);
  free(tuple2);
  
  ReleaseMemoryDispersion();
}

void initialize_arrays(int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    for (int k = 0; k < K; k++) {
        s_min_distance_per_cluster[k] = INFINITY;
        s_min_distance_tuple[k][0] = 0;
        s_min_distance_tuple[k][1] = 0;
    }
}

void fill_arrays(int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    initialize_arrays(s_min_distance_tuple, s_min_distance_per_cluster);
    // Function to fill the arrays based on the distance matrix and cluster assignments
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            if (Distances[i][j] < s_min_distance_per_cluster[partition[i]] && partition[i] == partition[j]) { // IMPORTANT: add "&& partition[i] = ..." to check that both data points are in the same cluster!
                s_min_distance_per_cluster[partition[i]] = Distances[i][j];
                s_min_distance_tuple[partition[i]][0] = i;
                s_min_distance_tuple[partition[i]][1] = j;
            }
        }
    }
}


// Function to evaluate the objective function
double evaluate_objective(double *s_min_distance_per_cluster) {
    double f = s_min_distance_per_cluster[0];
    for (int k = 1; k < K; k++) {
        f = fmin(f, s_min_distance_per_cluster[k]);
    }
    return f;
}

// convenience function that recalculates the min cluster distance of cluster k and changes the values in-place
void recalculate_cluster_distance(int k, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    s_min_distance_per_cluster[k] = INFINITY;
    for (int i = 0; i < N - 1; i++) {
        if (partition[i] != k) continue;
        for (int j = i + 1; j < N; j++) {
            if (Distances[i][j] < s_min_distance_per_cluster[k] && partition[j] == k) {
                s_min_distance_per_cluster[k] = Distances[i][j];
                s_min_distance_tuple[k][0] = i;
                s_min_distance_tuple[k][1] = j;
            }
        }
    }
}

// Function to add an element to a cluster
void adding(int new_ind, int cluster, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    for (int i = 0; i < N; i++) {
        if (partition[i] == cluster && i != new_ind) { // add "&& i != new_ind" to make sure to not include d[i][i] = 0!
            if (Distances[i][new_ind] < s_min_distance_per_cluster[cluster]) {
                s_min_distance_per_cluster[cluster] = Distances[i][new_ind];
                s_min_distance_tuple[cluster][0] = i;
                s_min_distance_tuple[cluster][1] = new_ind;
            }
        }
    }
    partition[new_ind] = cluster;
}

// Function to remove an element from its cluster
void removing(int removed_ind, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster)  {
    int g = partition[removed_ind];
    partition[removed_ind] = -1;  // Temporarily hide the index
    if (s_min_distance_tuple[g][0] == removed_ind || s_min_distance_tuple[g][1] == removed_ind) {
		// only requires update of objective function if remove_ind appears in min_distance_tuple
        s_min_distance_per_cluster[g] = INFINITY;
        for (int i = 0; i < N - 1; i++) {
            if (partition[i] == g) {
                for (int j = i + 1; j < N; j++) {
                    if (partition[j] == g) {
                        if (Distances[i][j] < s_min_distance_per_cluster[g]) {
                            s_min_distance_per_cluster[g] = Distances[i][j];
                            s_min_distance_tuple[g][0] = i;
                            s_min_distance_tuple[g][1] = j;
                        }
                    }
                }
            }
        }
    }
}

void swapping(int ind1, int ind2, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    int g1 = partition[ind1];
    int g2 = partition[ind2];

    // Temporarily hide ind1
    partition[ind1] = -1;

    // if either ind1/ind2 belong to their respective min_distance_tuple, ensure to recalculate the corresponding s_min_distance_per_cluster when adding ind1/ind2
   if (ind1 == s_min_distance_tuple[g1][0] || ind1 == s_min_distance_tuple[g1][1]) {
        recalculate_cluster_distance(g1,partition,s_min_distance_tuple,s_min_distance_per_cluster);
    }
    if (ind2 == s_min_distance_tuple[g2][0] || ind2 == s_min_distance_tuple[g2][1]) {
        partition[ind2] = -1;
        recalculate_cluster_distance(g2,partition,s_min_distance_tuple,s_min_distance_per_cluster);
        partition[ind2] = g2;
    }
    // the correct recalculation of s_min_distance_tuple happens in the adding-functions below

    // Add ind2 to the cluster of g1
    adding(ind2, g1, partition, s_min_distance_tuple, s_min_distance_per_cluster);

    // Add ind1 to the cluster of g2
    adding(ind1, g2, partition, s_min_distance_tuple, s_min_distance_per_cluster);
}


void SearchAlgorithmDisperion(void) {
    /* Algorithm 1: The main procedure of TPSDP. */

    int eta;
    int pickedSolution;

    //important! for windows and linux there is a differnt definition of this time
    //on windows its the wall time, on linux the CPU time
    clock_t start_time = clock();
    S_best.objective = -1;  //dispersion should be maximized
    
    // Initial population generation
    int i, j, k;
    for (i = 0; i < beta_max; i++) {

        /* Algorithm 2: initializes a solution S_D[i] */
        RandomInitialSol(S_D[i].s, S_D[i].SizeG);
        DoubleNeighborhoodLocalSearchDispersion(S_D[i].s, S_D[i].SizeG, &(S_D[i].objective));

        /* Update the best solution S_best if necessary */
        if (S_best.objective < S_D[i].objective) {
            for (j = 0; j < N; j++) S_best.s[j] = S_D[i].s[j];
            for (k = 0; k < K; k++) S_best.SizeG[k] = S_D[i].SizeG[k];
            S_best.objective = S_D[i].objective;
        }
    }
 
    for (int counter = 1; counter <= maxNumberIterations; counter++) {
        
        eta = (int)(theta * N / K);
        for (i = 0; i < beta_max; i++) {
            for (j = 0; j < N; j++) O_D[i].s[j] = S_D[i].s[j];
            for (k = 0; k < K; k++) O_D[i].SizeG[k] = S_D[i].SizeG[k];
            O_D[i].objective = S_D[i].objective;
        }
        // Strong Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            UndirectedPerturbation(eta, S_D[i].s, S_D[i].SizeG);
            DoubleNeighborhoodLocalSearchDispersion(S_D[i].s, S_D[i].SizeG, &S_D[i].objective);
        
            if (S_best.objective < S_D[i].objective) {
                for (j = 0; j < N; j++) S_best.s[j] = S_D[i].s[j];
                for (k = 0; k < K; k++) S_best.SizeG[k] = S_D[i].SizeG[k];
                S_best.objective = S_D[i].objective;
            }
        }

        // Crossover and Local Search
        if (beta_max > 1) {
            for (i = 0; i < beta_max; i++) {
                pickedSolution = random_int(beta_max);
                if (pickedSolution == i) {
                    pickedSolution = (pickedSolution + 1) % beta_max;
                }

                CrossoverDispersion(S_D[i].s, S_D[pickedSolution].s, O_D[i].s, O_D[i].SizeG);
                DoubleNeighborhoodLocalSearchDispersion(O_D[i].s, O_D[i].SizeG, &O_D[i].objective);
            }
            for (i = 0; i < beta_max; i++) {
                if (O_D[i].objective >= S_D[i].objective) {
                    for (j = 0; j < N; j++) S_D[i].s[j] = O_D[i].s[j];
                    for (k = 0; k < K; k++) S_D[i].SizeG[k] = O_D[i].SizeG[k];
                    S_D[i].objective = O_D[i].objective;
                } else if (LocalSearchCriterionCalculation(&O_D[i], &S_D[i]) > 1) {
                    for (j = 0; j < N; j++) S_D[i].s[j] = O_D[i].s[j];
                    for (k = 0; k < K; k++) S_D[i].SizeG[k] = O_D[i].SizeG[k];
                    S_D[i].objective = O_D[i].objective;
                }

                if (S_best.objective < S_D[i].objective) {
                    for (j = 0; j < N; j++) S_best.s[j] = S_D[i].s[j];
                    for (k = 0; k < K; k++) S_best.SizeG[k] = S_D[i].SizeG[k];
                    S_best.objective = S_D[i].objective;
                }
            }
        }

        // Direct Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            DirectPerturbationDispersion(eta_max, S_D[i].s, S_D[i].SizeG);
            DoubleNeighborhoodLocalSearchDispersion(S_D[i].s, S_D[i].SizeG, &S_D[i].objective);

            if (S_D[i].objective > S_best.objective) {
                for (j = 0; j < N; j++) S_best.s[j] = S_D[i].s[j];
                for (k = 0; k < K; k++) S_best.SizeG[k] = S_D[i].SizeG[k];
                S_best.objective = S_D[i].objective;
            }
        }

        // Linearly decrease population size
        qsort(S_D, beta_max, sizeof(Solution), CompareSolution);
        beta_max = (int)(beta_min - beta_max) * counter / maxNumberIterations + beta_max;
        theta = theta_max - (theta_max - theta_min) * counter / maxNumberIterations;
    }

    // Stop measuring time and calculate the elapsed time
    clock_t end_time = clock();
    double elapsed_time = (double) (end_time - start_time)/CLOCKS_PER_SEC;
    // Rprintf("The run time of the distance_clustering algortihm in seconds is: %f\n", elapsed_time);
}

void DoubleNeighborhoodLocalSearchDispersion(int s[], int SizeGroup[], double* objective) {

    int v, g, u;

    // Initialize the delta_f value and tuple arrays
    double delta_f = -99999.0;
    fill_arrays(s, min_distance_tuple, min_distance_per_cluster);

    int g1, g2;
    double old_f1, old_f2, delta_min;
    int imp;
    do { 
        imp = 0;  // Reset improvement flag

        // First loop: Move individual elements to improve partition
        for (v = 0; v < N; v++) {
            for (g = 0; g < K; g++) {
                // Check if moving `v` to `group` is valid and beneficial
                if ((s[v] != g) && (SizeGroup[s[v]] > LB[s[v]]) && (SizeGroup[g] < UB[g])) {
                    g1 = s[v];
                    old_f1 = min_distance_per_cluster[s[v]];
                    old_f2 = min_distance_per_cluster[g];
                    tuple1[0] = min_distance_tuple[g1][0]; tuple1[1] = min_distance_tuple[g1][1];
                    tuple2[0] = min_distance_tuple[g ][0]; tuple2[1] = min_distance_tuple[g ][1];
                    removing(v, s, min_distance_tuple, min_distance_per_cluster);
	                adding(v, g, s, min_distance_tuple, min_distance_per_cluster);

                     // accept changes only if the smaller of the two does not get down-graded (possibly decreasing the global minimum)
                    if (old_f1 < old_f2) delta_min = min_distance_per_cluster[g1] - old_f1;
                    else if (old_f1 > old_f2) delta_min = min_distance_per_cluster[g] - old_f2;
                    else delta_min = fmin(min_distance_per_cluster[g1] - old_f1,min_distance_per_cluster[g] - old_f2);

                    delta_f = min_distance_per_cluster[g] - old_f2 + min_distance_per_cluster[g1] - old_f1;
                    if (delta_f <= 0 || delta_min < 0) {
                        // revert changes
                        s[v] = g1;
                        min_distance_per_cluster[g1] = old_f1;
                        min_distance_per_cluster[g] = old_f2;
                        min_distance_tuple[g1][0] = tuple1[0]; min_distance_tuple[g1][1] = tuple1[1];
                        min_distance_tuple[g ][0] = tuple2[0]; min_distance_tuple[g ][1] = tuple2[1];
                    } else {
                        // update new group sizes 
                        SizeGroup[g]  += 1;
                        SizeGroup[g1] -= 1;
                        // maybe change this to a for loop 
                        // when set to 1, the group protenitally can run a long time
                        imp = 0;
                    }
                }
            }
        }

        // Second loop: Swap pairs of elements between groups
        for (v = 0; v < N-1; v++) {
            for (u = v + 1; u < N; u++) {
                // Only swap if nodes are in different groups
                if (s[v] != s[u]) {
                    g1 = s[v];
                    g2 = s[u];
                    old_f1 = min_distance_per_cluster[s[v]];
                    old_f2 = min_distance_per_cluster[s[u]];

                    tuple1[0] = min_distance_tuple[g1][0]; tuple1[1] = min_distance_tuple[g1][1]; // maybe quicker to just set tuple1 = min_distance_tuple[g1]? Does that work?
                    tuple2[0] = min_distance_tuple[g2][0]; tuple2[1] = min_distance_tuple[g2][1];
                    swapping(u, v, s, min_distance_tuple, min_distance_per_cluster); 

                    // accept changes only if the smaller of the two does not get down-graded (possibly decreasing the global minimum)
                    if (old_f1 < old_f2) delta_min = min_distance_per_cluster[g1] - old_f1;
                    else if (old_f1 > old_f2) delta_min = min_distance_per_cluster[g2] - old_f2;
                    else delta_min = fmin(min_distance_per_cluster[g1] - old_f1,min_distance_per_cluster[g2] - old_f2);

                    delta_f = min_distance_per_cluster[g2] - old_f2 + min_distance_per_cluster[g1] - old_f1; // IMPORTANT: typo: g2 instead of g
                    if (delta_f <= 0 || delta_min < 0) {
                        // revert changes
                        s[v] = g1;
                        s[u] = g2;
                        min_distance_per_cluster[g1] = old_f1;
                        min_distance_per_cluster[g2] = old_f2;
                        min_distance_tuple[g1][0] = tuple1[0]; min_distance_tuple[g1][1] = tuple1[1];
                        min_distance_tuple[g2][0] = tuple2[0]; min_distance_tuple[g2][1] = tuple2[1];
                    } else {
                        // maybe change this to a for loop 
                        // when set to 1, the group protenitally can run a long time
                        imp = 0;
                    }
                }
            }
        }
    } while (imp == 1);  // Continue until no improvement is made
    
    // Before writing objective to *objective, first update this value!
    f_objective = evaluate_objective(min_distance_per_cluster);
    
    *objective = f_objective;
}


void DirectPerturbationDispersion(int eta_max, int s[], int SizeGroup[]) {
    /* Algorithm 6: Directed Perturbation. 
	Iteratively refines partitions to balance group sizes and minimize costs */

    int k;
    int new_ind, ind1, ind2, selectedElement, minElement;
    double objective_k, minDeltaValue, maxDeltaValue;

    fill_arrays(s, min_distance_tuple, min_distance_per_cluster);

    // Main loop for perturbation iterations
    int number;
    for (int L = 0; L < eta_max; L++) {

        // Reset tracking variables
        number = 0;
        for (int i = 0; i < K; i++) {
            UnderLB[i] = 0;
            Rd[i] = -1;
        }

        for (k = 0; k < K; k++) {
            // Pseudo group 8 to 13
            /* Remove the element x, y from the tuple of cluster k which have the lowest dipserion.
             Calculate the minimal dispersion value after removing x and y in cluster k. 
             The element which leads to the lower minimum dispersion of the cluster k after its removal
             will be removed from the cluster k. */
          
            // store values before removing() in temporary values. This speeds-up the adding-process. In short:
            ind1 = min_distance_tuple[k][0];
            ind2 = min_distance_tuple[k][1];
            objective_k = min_distance_per_cluster[k];
            removing(ind1, s, min_distance_tuple, min_distance_per_cluster);
            minDeltaValue = min_distance_per_cluster[k];
            // restore changes simply from ind1 and objective_k:
            min_distance_tuple[k][0] = ind1;
            min_distance_tuple[k][1] = ind2;
            min_distance_per_cluster[k] = objective_k;
            s[ind1] = k; // restore class membership
            removing(ind2, s, min_distance_tuple, min_distance_per_cluster);
             if (min_distance_per_cluster[k] < minDeltaValue) {
                minDeltaValue = min_distance_per_cluster[k];
                minElement = ind2;
            } else {
                min_distance_tuple[k][0] = ind1;
                min_distance_tuple[k][1] = ind2;
                min_distance_per_cluster[k] = objective_k;
                s[ind2] = k; // restore class membership
                removing(ind1, s, min_distance_tuple, min_distance_per_cluster);
                minElement = ind1;
            }

            // Record the minimum element for removal
            Rd[k] = minElement;
            SizeGroup[k] -= 1;

            // If the group size falls below the lower bound, mark it
            if (SizeGroup[k] < LB[k]) {
                UnderLB[k] = 1;
                number += 1;
            }
        }

        // line 15 of pseudo group will be removed, since it is not necessary for dispersion     
		
        // Handle groups that are under the lower bound (LB)
        int selectedGroup;
        int nn = 0;
        while (nn < number) {
            // pseudo code line 17 - 29
            k = random_int(K);

            // Find the element with the highest average connection to the group
            while (UnderLB[k] == 0) {
                k = (k + 1) % K;
            }
            
            // IMPORTANT: add the next five lines here
            // IF cluster k has less than 2 elements its min_distance is set to INFINITY.
            // In this case, we can safely set it to 0
            objective_k = min_distance_per_cluster[k];
            bool is_infinite = objective_k == INFINITY;
            if (is_infinite) objective_k = 0.0;

            // store values before adding() in temporary values. This speeds-up the removing-process. In short:
            maxDeltaValue = -INFINITY; // min_delta cannot become positive by adding a new data point to a cluster for dispersion
            for (int i = 0; i < K; i++) {
                new_ind = Rd[i];
                if (new_ind > -1) {
                    ind1 = min_distance_tuple[k][0];
                    ind2 = min_distance_tuple[k][1];
                    adding(new_ind, k, s, min_distance_tuple, min_distance_per_cluster);
                    if (min_distance_per_cluster[k] - objective_k > maxDeltaValue) {
                        maxDeltaValue = min_distance_per_cluster[k] - objective_k;
                        selectedElement = new_ind;
                        selectedGroup = i;
                    }
                    // restore changes simply from ind1 and objective_k:
                    min_distance_tuple[k][0] = ind1;
                    min_distance_tuple[k][1] = ind2;
                    // replace "min_distance_per_cluster[k] = objective_k;" by the if-else statement, see other IMPORTANT above
                    if (is_infinite) min_distance_per_cluster[k] = INFINITY;
                    else min_distance_per_cluster[k] = objective_k;
                    s[new_ind] = -1; // remove new_ind from class
                }
            }

            adding(selectedElement, k, s, min_distance_tuple, min_distance_per_cluster);  
            
            UnderLB[k] = 0;
            Rd[selectedGroup] = -1;
            nn++;
        }

        // Handle remaining elements for groups that are above LB
        nn = 0;
        while (nn < K - number) {
            // pseuo group line 32 - 43
            selectedGroup = random_int(K);
            while (Rd[selectedGroup] == -1) {
                selectedGroup = (selectedGroup + 1) % K;
            }

            maxDeltaValue = -INFINITY;
            new_ind = Rd[selectedGroup];
            Rd[selectedGroup] = -1;
            for (k = 0; k < K; k++) {
                // check that upper-bound tmpUB[k] is not yet reached
                    if (SizeGroup[k] < UB[k]) {
                    ind1 = min_distance_tuple[k][0];
                    ind2 = min_distance_tuple[k][1];
                    objective_k = min_distance_per_cluster[k];
                    adding(new_ind, k, s, min_distance_tuple, min_distance_per_cluster);
                    if (min_distance_per_cluster[k] - objective_k > maxDeltaValue) {
                        minDeltaValue = min_distance_per_cluster[k] - objective_k;
                        selectedGroup = k;
                    }
                    // restore changes simply from ind1 and objective_k:
                    min_distance_tuple[k][0] = ind1;
                    min_distance_tuple[k][1] = ind2;
                    min_distance_per_cluster[k] = objective_k;
                    s[new_ind] = -1; // remove new_ind from class
                }
            }
            adding(new_ind, selectedGroup, s, min_distance_tuple, min_distance_per_cluster); 

            nn += 1;
        }
    }
}

void CrossoverDispersion(int partition1[], int partition2[], int solutionChild[], int scSizeGroup[]) {
    /* Algorithm 5: combines partitions in a way that maintains group constraints */

    int i, j, selectedGroup;
    double maxGroupDispersion;
    int elementCount, groupCount;
    int targetGroup = -1;
    int processedCount;
    int selectedElement;
    int totalLowerBound, totalBelowLowerBound;

    // Initialize s1 and p2s with partition1
    for (i = 0; i < N; i++) {
        s1[i] = partition1[i];
        s2[i] = partition2[i];
    }

    fill_arrays(s1, min_distance_tuple_s1, min_distance_per_cluster_s1);
    fill_arrays(s2, min_distance_tuple_s2, min_distance_per_cluster_s2);
    for (int k = 0; k < K; k++) {
        min_distance_per_cluster_s1[k] = min_distance_per_cluster[k];  //groupDiversity
        min_distance_per_cluster_s2[k] = min_distance_per_cluster[k];
        min_distance_tuple_s1[k][0] = min_distance_tuple[k][0];
        min_distance_tuple_s2[k][1] = min_distance_tuple[k][0];
    }
    
    // Initialize arrays
    for (i = 0; i < N; i++) {
        vectorElement[i] = i;
        solutionChild[i] = -1; // ?
    }
    for (i = 0; i < K; i++) {
        LBGroup[i] = 0;
        UBGroup[i] = 0;
        BigThanLB[i] = 0;
        groupElement[i] = i;
        tmpUB[i] = UB[i];
        scSizeGroup[i] = 0;
    }

    // Main crossover process
    for (i = 0; i < K; i++) {
        if (uniform_rnd_number() < 0.5) {
            // Process partition1
            //find group with highest dispersion
            maxGroupDispersion = -1;
            for (j = 0; j < K; j++) {
                if (min_distance_per_cluster_s1[j] > maxGroupDispersion) {
                    maxGroupDispersion = min_distance_per_cluster_s1[j];
                    selectedGroup = j;
                }
            }

            // select elements to move
            elementCount = 0;
            for (j = 0; j < N; j++ ) {
                if (s1[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

           // choose groups who have enough space to hole seleccted groups
            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (tmpUB[j] != -1 && tmpUB[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { 
                // 16  from Pseudogroup Algotihm 5 
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (tmpUB[j] != -1 && elementCount - tmpUB[j] < minDiff) {
                        minDiff = elementCount - tmpUB[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    // 17  from Pseudogroup Algotihm 5 
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    solutionChild[SelectEle[selectedElement]] = targetGroup;
                    tmpEle[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];  // 13  from Pseudogroup Algotihm 5 
                for (j = 0; j < elementCount; j++) {
                     // 14  from Pseudogroup Algotihm 5 
                    solutionChild[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    tmpEle[j] = SelectEle[j];
                }
            }
        } else {
            // Process partition2 
            maxGroupDispersion = -1;
            for (j = 0; j < K; j++) {
                if (min_distance_per_cluster_s2[j] > maxGroupDispersion) {
                    maxGroupDispersion = min_distance_per_cluster_s2[j];
                    selectedGroup = j;
                }
            }

            elementCount = 0;
            for (j = 0; j < N; j++) {
                if (s2[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (tmpUB[j] != -1 && tmpUB[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { // No valid group found
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (tmpUB[j] != -1 && elementCount - tmpUB[j] < minDiff) {
                        minDiff = elementCount - tmpUB[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    solutionChild[SelectEle[selectedElement]] = targetGroup;
                    tmpEle[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];
                for (j = 0; j < elementCount; j++) {
                    solutionChild[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    tmpEle[j] = SelectEle[j];
                }
            }
        }

        // Update group dispersion
        for (j = 0; j < elementCount; j++) {
            removing(tmpEle[j], s1, min_distance_tuple_s1, min_distance_per_cluster_s1);
            removing(tmpEle[j], s2, min_distance_tuple_s2, min_distance_per_cluster_s2);
            s1[tmpEle[j]] = -1;
            s2[tmpEle[j]] = -1;
        }

        // group should not be available anymore
        min_distance_per_cluster_s1[targetGroup] = -1;
        min_distance_per_cluster_s2[targetGroup] = -1;
        tmpUB[targetGroup] = -1;
        scSizeGroup[targetGroup] = elementCount;
    }

    // Adjust assignments to maintain group size constraints
    processedCount = 0;
    totalLowerBound = 0;
    totalBelowLowerBound = 0;
    for (i = 0; i < K; i++) {
        totalLowerBound += LB[i];
        if (scSizeGroup[i] < LB[i]) {
            processedCount += scSizeGroup[i];
            totalBelowLowerBound += scSizeGroup[i];
            LBGroup[i] = 1;
        } else {
            processedCount += LB[i];
        }
        if (scSizeGroup[i] > LB[i]) {
            BigThanLB[i] = 1;
        }
    }

    // Assign unprocessed elements
    for (i = 0; i < N; i++) {
        if (vectorElement[i] != -1) {
            processedCount++;
        }
    }

    while (processedCount < totalLowerBound) {
        targetGroup = random_int(K);
        while (BigThanLB[targetGroup] == 0) {
            targetGroup = (targetGroup + 1) % K;
        } 

        elementCount = 0;
        for (j = 0; j < N; j++) {
            if (solutionChild[j] == targetGroup) {
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = random_int(elementCount);
        solutionChild[SelectEle[selectedElement]] = -1;
        vectorElement[SelectEle[selectedElement]] = SelectEle[selectedElement];
        scSizeGroup[targetGroup]--;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            BigThanLB[targetGroup] = 0;
        }
        processedCount++;
    }

    int sumLB = 0;
    for (i = 0; i < K; i++) {
        if (LBGroup[i] == 1) {
            sumLB += LB[i];
        }
    }

    while (totalBelowLowerBound < sumLB) {
        targetGroup = random_int(K);
        while (LBGroup[targetGroup] == 0) {
            targetGroup = (targetGroup + 1) % K;
        }

        elementCount = 0;
        for (i = 0; i < N; i++) {
            if (vectorElement[i] != -1) {
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = random_int(elementCount);
        solutionChild[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            LBGroup[targetGroup] = 0;
        }
        totalBelowLowerBound++;
    }

    int totalSize = 0;
    for (i = 0; i < K; i++) {
        totalSize += scSizeGroup[i];
        if (scSizeGroup[i] < UB[i]) {
            UBGroup[i] = 1;
        }
    }

    while (totalSize < N) {
        targetGroup = random_int(K);
        while (UBGroup[targetGroup] == 0) {
            targetGroup = (targetGroup + 1) % K;
        }

        elementCount = 0;
        for (i = 0; i < N; i++) {
            if (vectorElement[i] != -1) {
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = random_int(elementCount);
        solutionChild[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == UB[targetGroup]) {
            UBGroup[targetGroup] = 0;
        }
        totalSize++;
    }
}

void AssignMemoryDispersion(void) {
    /*  Allocates memory dynamically for various arrays and matrices necessary 
	for the algorithm's execution. This includes structures for population management, 
	distance matrices, diversity measures, and neighborhood exploration.
	*/
    
    S_D = (Solution*)malloc(beta_max * sizeof(Solution));
    O_D = (Solution*)malloc(beta_max * sizeof(Solution));
    
    int i; 
    for (i = 0; i < beta_max; i++) {
        S_D[i].s = (int*)malloc(N * sizeof(int));
        O_D[i].s = (int*)malloc(N * sizeof(int));
        S_D[i].SizeG = (int*)malloc(K * sizeof(int));
        O_D[i].SizeG = (int*)malloc(K * sizeof(int));
    }
    
    S_best.s = (int*)malloc(N * sizeof(int));
    S_best.SizeG = (int*)malloc(K * sizeof(int));
    
    Rd = (int*)malloc(K * sizeof(int));
    for (i = 0; i < K; i++) Rd[i] = 0;
    UnderLB = (int*)malloc(K * sizeof(int));
    
    tmpUB = (int*)malloc(K * sizeof(int));
    LBGroup = (int*)malloc(K * sizeof(int));
    UBGroup = (int*)malloc(K * sizeof(int));
    BigThanLB = (int*)malloc(K * sizeof(int));
    vectorElement = (int*)malloc(N * sizeof(int));
    groupElement = (int*)malloc(K * sizeof(int));
    SelectEle = (int*)malloc(N * sizeof(int));
    SelectGroup = (int*)malloc(K * sizeof(int));
    tmpEle = (int*)malloc(N * sizeof(int));
    s1 = (int*)malloc(N * sizeof(int));
    s2 = (int*)malloc(N * sizeof(int));
}

void ReleaseMemoryDispersion(void) {
    /* responsible for reading the input file, 
    initializing matrices, and setting constraints on group sizes. */ 

    free(S_best.s); S_best.s = NULL;
    free(S_best.SizeG); S_best.SizeG = NULL;
    
   // Rprintf("relesee dispersion Until now runs trhough.");
    // IMPORTANT: releasing S_D and O_D like for TPSDP leads to error!
    // int i;
    //for (i = 0; i < beta_max; i++) {
      //  free(S_D[i].s); S_D[i].s = NULL;
       // free(S_D[i].SizeG); S_D[i].SizeG = NULL;
       // free(O_D[i].s); O_D[i].s = NULL;
        //free(O_D[i].SizeG); O_D[i].SizeG = NULL;
   // }
    free(O_D); O_D = NULL;
    free(S_D); S_D = NULL; 
   // Rprintf("relesee dispersion Until now runs trhough 2.");

       
    free(LB); LB = NULL;
    free(UB); UB = NULL;

    free(Rd); Rd = NULL;
    free(UnderLB); UnderLB = NULL;
    free(tmpUB); tmpUB = NULL;
    free(LBGroup); LBGroup = NULL;
    free(UBGroup); UBGroup = NULL;
    free(BigThanLB); BigThanLB = NULL;
    free(vectorElement); vectorElement = NULL;
    free(groupElement); groupElement = NULL;
    free(SelectEle); SelectEle = NULL;
    free(SelectGroup); SelectGroup = NULL;
    free(tmpEle); tmpEle = NULL;
    free(s1); s1 = NULL;
    free(s2); s2 = NULL;
}
