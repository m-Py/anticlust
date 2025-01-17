#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "three-phase-header.h" 

int beta_max;
int N, K;  // node number and group number
double** Distances;   // distance matrix between elements
double** DistancesT;
int* LB; // Lower bound of the number of elements in the group i 
int* UB; // Upper bound of the number of elements in the group i
double theta, theta_max, theta_min; 
double alpha;
int beta_min; 
int eta_max;

// Matrix M
double **Delta_Matrix;  // incremental matrix 
double **Delta_Matrix_p1;
double **Delta_Matrix_p2;

// main processure
int maxNumberIterations;
double start_time, Time_limit;
Solution S_best; //best solution
Solution *S; //S_i
Solution *O; //O_i in crossover

//double neighboorhood local search
double f_objective;

// for crossover
int *vectorElement;
int *groupElement;
int *SelectGroup;
int *SelectEle;
int *tmpEle;
int *s1;
int *s2;
double *groupDiversity_s1;
double *groupDiversity_s2;
int *LBGroup, *UBGroup, *tmpUB, *BigThanLB;

//directed pertubation
int* Rd, *UnderLB; //Rd=R
double** Avg;

void DirectPerturbationDiversity(int eta_max, int partition[], int SizeG[]);
void ClearDeltaMatrix(void);
void BuildDeltaMatrix(int partition[]);
void OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup);
void BuildGroupDiversityForCrossover(int partition[], double groupDiversity[]);
void AssignMemoryDiversity(void);
void ReleaseMemoryDiversity(void);
void SearchAlgorithmDiversity(void);
void DoubleNeighborhoodLocalSearchDiversity(int partition[], int SizeGroup[], double* objective);
void CrossoverDiversity(int partition1[], int partition2[], int score[], int scSizeGroup[]);

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
void three_phase_search_dynamic_population_size(
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
  
  // Fill Distances and DistancesT with values from input
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Distances[i][j] = distances[i * N + j];
      DistancesT[i][j] = 2 * distances[i * N + j];
    }
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
  
  AssignMemoryDiversity();
  if (*mem_error == 1) {
    return;
  }
    
  SearchAlgorithmDiversity();
  
  //save S_best -> solution with result
  for (int i = 0; i < N; i++){
    result[i] = S_best.s[i];
  }
  *score = S_best.objective;
  *elapsed_time = Time_limit;
  
  // Remember to free the allocated memory after use
  for (int i = 0; i < N; i++) {
    free(Distances[i]); Distances[i] = NULL;
    free(DistancesT[i]); DistancesT[i] = NULL;
  }
  free(Distances); Distances = NULL;
  free(DistancesT); DistancesT = NULL;
  free(LB); LB = NULL;
  free(UB); UB = NULL;
  
  ReleaseMemoryDiversity();
}


void SearchAlgorithmDiversity(void) {
    /* Algorithm 1: The main procedure of TPSDP. */

    int eta;
    int pickedSolution;

    //important! for windows and linux there is a differnt definition of this time
    //on windows its the wall time, on linux the CPU time
    clock_t start_time = clock();
    S_best.objective = -INFINITY;
    
    // Initial population generation
    int i, j, k;
    for (i = 0; i < beta_max; i++) {
        /* Algorithm 2: initializes a solution S_D[i] */
        RandomInitialSol(S[i].s, S[i].SizeG);
        DoubleNeighborhoodLocalSearchDiversity(S[i].s, S[i].SizeG, &(S[i].objective));

        /* Update the best solution S_best if necessary */
        if (S_best.objective < S[i].objective) {
            for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
            for (k = 0; k < K; k++) S_best.SizeG[k] = S[i].SizeG[k];
            S_best.objective = S[i].objective;
        }
    }
 
    for (int counter = 1; counter <= maxNumberIterations; counter++) {
        
        eta = (int)(theta * N / K);
        for (i = 0; i < beta_max; i++) {
            for (j = 0; j < N; j++) O[i].s[j] = S[i].s[j];
            for (k = 0; k < K; k++) O[i].SizeG[k] = S[i].SizeG[k];
            O[i].objective = S[i].objective;
        }
        // Strong Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            UndirectedPerturbation(eta, S[i].s, S[i].SizeG);
            DoubleNeighborhoodLocalSearchDiversity(S[i].s, S[i].SizeG, &S[i].objective);

            if (S_best.objective < S[i].objective) {
                for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                for (k = 0; k < K; k++) S_best.SizeG[k] = S[i].SizeG[k];
                S_best.objective = S[i].objective;
            }
        }

        // Crossover and Local Search
        if (beta_max > 1) {
            for (i = 0; i < beta_max; i++) {
                pickedSolution = random_int(beta_max);
                if (pickedSolution == i) {
                    pickedSolution = (pickedSolution + 1) % beta_max;
                }

                CrossoverDiversity(S[i].s, S[pickedSolution].s, O[i].s, O[i].SizeG);
                DoubleNeighborhoodLocalSearchDiversity(O[i].s, O[i].SizeG, &O[i].objective);
            }
            for (i = 0; i < beta_max; i++) {
                if (O[i].objective >= S[i].objective) {
                    for (j = 0; j < N; j++) S[i].s[j] = O[i].s[j];
                    for (k = 0; k < K; k++) S[i].SizeG[k] = O[i].SizeG[k];
                    S[i].objective = O[i].objective;
                } else if (LocalSearchCriterionCalculation(&O[i], &S[i]) > 1) {
                    for (j = 0; j < N; j++) S[i].s[j] = O[i].s[j];
                    for (k = 0; k < K; k++) S[i].SizeG[k] = O[i].SizeG[k];
                    S[i].objective = O[i].objective;
                }
                
                if (S_best.objective < S[i].objective) {
                    for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                    for (k = 0; k < K; k++) S_best.SizeG[k] = S[i].SizeG[k];
                    S_best.objective = S[i].objective;
                }
            }
        }

        // Direct Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            DirectPerturbationDiversity(eta_max, S[i].s, S[i].SizeG);
            DoubleNeighborhoodLocalSearchDiversity(S[i].s, S[i].SizeG, &S[i].objective);

            if (S[i].objective > S_best.objective) {
                for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                for (k = 0; k < K; k++) S_best.SizeG[k] = S[i].SizeG[k];
                S_best.objective = S[i].objective;
            }
        }

        // Linearly decrease population size
        qsort(S, beta_max, sizeof(Solution), CompareSolution);
        beta_max = (int)(beta_min - beta_max) * counter / maxNumberIterations + beta_max;
        theta = theta_max - (theta_max - theta_min) * counter / maxNumberIterations;
    }

    // Stop measuring time and calculate the elapsed time
    clock_t end_time = clock();
    double elapsed_time = (double) (end_time - start_time)/CLOCKS_PER_SEC;
    Time_limit = elapsed_time;
}

int CompareSolution(const void *first, const void *second) {
    /* Compares two solutions based on their objective */
    Solution *solution1 = (Solution *)first;
    Solution *solution2 = (Solution *)second;

    if (solution1->objective < solution2->objective) {
        return 1;  // solution2 has a greater objective
    } else if (solution1->objective > solution2->objective) {
        return -1; // solution1 has a greater objective
    } else {
        return 0;  // Both costs are equal
    }
}

// Function to swap two elements
void swap_elements(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}


void fisher_yates_shuffle(int arr[], int n) {
    // Fisher-Yates shuffle to generate a random permutation 
    for (int i = n - 1; i > 0; i--) {
        // Generate a random index j where 0 <= j <= i
        int j = random_int(i + 1);

        // Swap arr[i] and arr[j]
        swap_elements(&arr[i], &arr[j]);
    }
}

void RandomInitialSol(int s[], int SizeG[]) {
	/* Algorithm 2: initializes a random solution that respects the group size constraints */

    // Allocates memory
	int* isAssigned = (int *)malloc(N * sizeof(int));  // Tracks if an element is assigned
	int* groupSize =(int *)malloc(K * sizeof(int)); // Stores the size of each group
    int* permutedIndexList =(int *)malloc(N * sizeof(int)); // Stores the size of each group
    int* permutedGroupList =(int *)malloc(K * sizeof(int)); // Stores the size of each group
	
    int i;
	for (i = 0; i < K; i++) groupSize[i] = 0;
	for (i = 0; i < N; i++) isAssigned[i] = 0;
    for (i = 0; i < N; i++) permutedIndexList[i] = i;
    for (i = 0; i < K; i++) permutedGroupList[i] = i;

    fisher_yates_shuffle(permutedIndexList, N);
    fisher_yates_shuffle(permutedGroupList, K);

    // Calculate the total number of elements that need to satisfy the lower bounds
    int total_assigned = 0;
    int total_LB = 0;
    for (int i = 0; i < K; i++) {
        total_LB += LB[i];
    }

    // First phase: Assign elements to satisfy lower bound constraints (LB)
    int selected_element = 0;
    while (total_assigned < total_LB) {

        for (int group = 0; group < K; group++) {
            if (groupSize[group] < LB[group]) {
                s[permutedIndexList[selected_element]] = group;
                isAssigned[permutedIndexList[selected_element]] = 1;
                groupSize[group]++;
                total_assigned++;
                break;  // Move to the next element once assigned
            }
        }
        selected_element++;
    }

	// Second phase: Assign remaining elements, respecting the upper bound (UB)
    while (total_assigned < N) {
        for (int group = 0; group < K; group++) {
            if (groupSize[permutedGroupList[group]] < UB[permutedGroupList[group]]) {
                s[permutedIndexList[selected_element]] = permutedGroupList[group];
                isAssigned[permutedIndexList[selected_element]] = 1;
                groupSize[permutedGroupList[group]]++;
                total_assigned++;
                break;
            }
        }
        fisher_yates_shuffle(permutedGroupList, K);
        selected_element++;
    }

    // Copy the final group sizes into the output array SizeG
   	for (i = 0; i < K; i++)  SizeG[i] = groupSize[i];

    // Free allocated memory
	free(groupSize); groupSize = NULL;
	free(isAssigned); isAssigned = NULL;
    free(permutedIndexList); permutedIndexList = NULL; 
    free(permutedGroupList); permutedGroupList = NULL;
}


void DoubleNeighborhoodLocalSearchDiversity(int s[], int SizeGroup[], double* objective) {
    const double DELTA_THRESHOLD = 0.0001;  // Define a constant for comparison threshold
    int v, g, u;
    int oldGroup, oldGroup1, t;

    // Build the delta_f matrix for objective changes
    BuildDeltaMatrix(s);

    // Initialize the delta_f value
    double delta_f = -99999.0;

    int imp;
    do {
        imp = 0;  // Reset improvement flag

        // First loop: Move individual elements to improve partition
        for (v = 0; v < N; v++) {
            for (g = 0; g < K; g++) {
                // Check if moving `v` to `group` is valid and beneficial
                if ((s[v] != g) && (SizeGroup[s[v]] > LB[s[v]]) && (SizeGroup[g] < UB[g])) {
                    delta_f = Delta_Matrix[v][g] - Delta_Matrix[v][s[v]];

                    if (delta_f > DELTA_THRESHOLD) {
                        oldGroup = s[v];

                        // Update delta_f matrix for the move
                        OneMoveUpdateDeltaMatrix(v, oldGroup, g);

                        // Update group sizes
                        SizeGroup[oldGroup] -= 1;
                        SizeGroup[g] += 1;

                        // Assign v to new group
                        s[v] = g;

                        // Update total objective
                        f_objective += delta_f;

                        // Mark as improved
                        imp = 1;
                    }
                }
            }
        }

        // Second loop: Swap pairs of elements between groups
        for (v = 0; v < N; v++) {
            for (u = v + 1; u < N; u++) {
                // Only swap if nodes are in different groups
                if (s[v] != s[u]) {
                    delta_f = (Delta_Matrix[v][s[u]] - Delta_Matrix[v][s[v]])
                          + (Delta_Matrix[u][s[v]] - Delta_Matrix[u][s[u]])
                          - DistancesT[v][u];

                    if (delta_f > DELTA_THRESHOLD) {
                        oldGroup = s[v];
                        oldGroup1 = s[u];

                        // Update delta_f matrix M for the swap
                        OneMoveUpdateDeltaMatrix(v, oldGroup, oldGroup1);
                        OneMoveUpdateDeltaMatrix(u, oldGroup1, oldGroup);

                        // Swap the two nodes between groups
                        t = s[v];
                        s[v] = s[u];
                        s[u] = t;

                        // Update total objective
                        f_objective += delta_f;

                        // Mark as improved
                        imp = 1;
                    }
                }
            }
        }
    } while (imp == 1);  // Continue until no improvement is made

    // Update the partition array with the final assignments
    *objective = f_objective;
}

void UndirectedPerturbation(int theta, int s[], int SizeGroup[]) {
    /* Algorithm 4: Undirected Perturbation. Applies a strong perturbation to the partition */

    int perturb_type;
    int v, g, x, y;
    int oldGroup, swap;

    int count = 0;
    int NumberNeighbors = N * (N - 1) / 2 + N * K;
     while (count < theta) {
        perturb_type = random_int(NumberNeighbors);

        if (perturb_type  < N * K) {  // Type 1: Random (element, group) perturbation
            v = random_int(N); // Randomly choose an element v
            g = random_int(K); // Randomly choose a group g

             if (s[v] != g && SizeGroup[s[v]] > LB[s[v]] && SizeGroup[g] < UB[g]) {
                oldGroup = s[v];
                SizeGroup[oldGroup]--;
                SizeGroup[g]++;
                s[v] = g;
                count++;
            }
        } 
        else { // Type 2: Random (element x, element y) perturbation
            x = random_int(N); // Randomly choose element x
            y = random_int(N); // Randomly choose element y

            // Apply perturbation if elements are in different groups
            if (s[x] != s[y] && x != y) {
                swap = s[x];
                s[x] = s[y];
                s[y] = swap;
                count++;
            }
        }
    }
}

void DirectPerturbationDiversity(int eta_max, int s[], int SizeGroup[]) {
    /* Algorithm 6: Directed Perturbation. 
	Iteratively refines partitions to balance group sizes and minimize costs */

    int i, j, k;
    int minDeltaValue, minElement;

    BuildDeltaMatrix(s);

    int number;
    // Main loop for perturbation iterations
    for (int L = 0; L < eta_max; L++) {

        // Reset tracking variables
        number = 0;
        for (i = 0; i < K; i++) {
            UnderLB[i] = 0;
            Rd[i] = -1;
            for (j = 0; j < K; j++) {
                Avg[i][j] = 0.0;
            }
        }

        // Find the minimum scoring element for each group
        for (k = 0; k < K; k++) {
            minDeltaValue = 99999999;
            minElement = 0;
            for (i = 0; i < N; i++) {
                if (s[i] == k) {
                    if (Delta_Matrix[i][k] < minDeltaValue) {
                        minDeltaValue = Delta_Matrix[i][k];
                        minElement = i;
                    }
                }
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

        // Rebuild the Delta matrix and average connections after removal
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                Delta_Matrix[Rd[i]][s[Rd[j]]] = Delta_Matrix[Rd[i]][s[Rd[j]]] - Distances[Rd[i]][Rd[j]];
                Avg[s[Rd[i]]][s[Rd[j]]] = Delta_Matrix[Rd[i]][s[Rd[j]]] / SizeGroup[s[Rd[j]]];
            }
        }        
		
        // Handle groups that are under the lower bound (LB)
        int selectedGroup;
        int maxAvgCon;
        int nn = 0;
        int i;
        while (nn < number) {
            maxAvgCon = -9999;
            i = random_int(K);

            // Find the element with the highest average connection to the group
            while (UnderLB[i] == 0) {
                i = (i + 1) % K;
            }
            for (j = 0; j < K; j++) {
                if (Avg[j][i] > maxAvgCon && Rd[j]!=-1) {
                    maxAvgCon = Avg[j][i];
                    selectedGroup = j;
                }
            }

            // Move the selected element to the group
            SizeGroup[i] += 1;
            for (k = 0; k < K; k++) {
                if (Rd[k] != -1) {
                    Delta_Matrix[Rd[k]][i] += Distances[Rd[k]][Rd[selectedGroup]];
                    Avg[s[Rd[k]]][i] = Delta_Matrix[Rd[k]][i] / SizeGroup[i];
                }
            }

            // Clear old connections for the moved element and finalize the move
            for (k = 0; k < K; k++) {
                Avg[s[Rd[selectedGroup]]][k] = 0.0;
            }
            s[Rd[selectedGroup]] = i;
            UnderLB[i] = 0;
            Rd[selectedGroup] = -1;
            nn++;
        }

        // Handle remaining elements for groups that are above LB
        int groupWithMaxAvgCon;
        nn = 0;
        while (nn < K - number) {
            selectedGroup = random_int(K);
            while (Rd[selectedGroup] == -1) {
                selectedGroup = (selectedGroup + 1) % K;
            }
            maxAvgCon = -9999;
            for (j = 0; j < K; j++) {
                if (Avg[selectedGroup][j] > maxAvgCon) {
                    maxAvgCon = Avg[selectedGroup][j];
                    groupWithMaxAvgCon = j;
                }
            }
            // Move the selected element to the group
            if (SizeGroup[groupWithMaxAvgCon] < UB[groupWithMaxAvgCon]) {
                SizeGroup[groupWithMaxAvgCon] += 1;
                for (k = 0; k < K; k++) {
                    if (Rd[k] != -1) {
                        Delta_Matrix[Rd[k]][groupWithMaxAvgCon] += Distances[Rd[k]][Rd[selectedGroup]];
                        Avg[s[Rd[k]]][groupWithMaxAvgCon] = Delta_Matrix[Rd[k]][groupWithMaxAvgCon] / SizeGroup[groupWithMaxAvgCon];
                    }
                }
                for (k = 0; k < K; k++) {
                	Avg[s[Rd[selectedGroup]]][k] = 0.0;
                }
                s[Rd[selectedGroup]] = groupWithMaxAvgCon;
				Rd[selectedGroup] = -1;
				nn += 1;
			}
			else {
				for (k = 0; k < K; k++) {
					Avg[k][groupWithMaxAvgCon] = 0.0;
				}
            }
        }
        BuildDeltaMatrix(s);
    }
}

// Function to process a partition 
void process_partition(double* groupDiversity, int* partition,  int* tmpUB, int* childSolution,
     int* vectorElement, int K, int N, int *element_count, int *target_group) {
    int i, selectedGroup, processedCount, selectedElement;

    int elementCount = *element_count; 
    int targetGroup = *target_group;
    int maxGroupDiversity = -1;
    for (i = 0; i < K; i++) {
        if (groupDiversity_s1[i] > maxGroupDiversity) {
            maxGroupDiversity = groupDiversity_s1[i];
            selectedGroup = i;
        }   
    }

    for (i = 0; i < N; i++) {
        if (s1[i] == selectedGroup) {
            SelectEle[elementCount++] = i;
        }
    }

    int groupCount = 0;
    for (i = 0; i < K; i++) {
        if (tmpUB[i] != -1 && tmpUB[i] >= elementCount) {
            SelectGroup[groupCount++] = i;
        }
    }

    if (groupCount == 0) { // No valid group found
        int minDiff = 999999;
        for (i = 0; i < K; i++) {
            if (tmpUB[i] != -1 && elementCount - tmpUB[i] < minDiff) {
                minDiff = elementCount - tmpUB[i];
                targetGroup = i;
            }
        }

        processedCount = 0;
        while (processedCount < elementCount - minDiff) {
            selectedElement = random_int(elementCount);
            while (SelectEle[selectedElement] == -1) {
                selectedElement = (selectedElement + 1) % elementCount;
            }

            childSolution[SelectEle[selectedElement]] = targetGroup;
            tmpEle[processedCount++] = SelectEle[selectedElement];
            vectorElement[SelectEle[selectedElement]] = -1;
            SelectEle[selectedElement] = -1;
        }
        elementCount = processedCount;
    } else {
        targetGroup = SelectGroup[random_int(groupCount)];
        for (i = 0; i < elementCount; i++) {
            childSolution[SelectEle[i]] = targetGroup;
            vectorElement[SelectEle[i]] = -1;
            tmpEle[i] = SelectEle[i];
        }
    }
    *element_count = elementCount;
    *target_group = targetGroup;
}

void CrossoverDiversity(int partition1[], int partition2[], int childSolution[], int scSizeGroup[]) {
    /* Algorithm 5: combines partitions in a way that maintains group constraints */

    int i, j;
    int elementCount, processedCount, selectedElement;
    int totalLowerBound, totalBelowLowerBound;

    // Initialize arrays
    for (i = 0; i < N; i++) {
        vectorElement[i] = i;
        childSolution[i] = -1;
    }
    for (i = 0; i < K; i++) {
        LBGroup[i] = 0;
        UBGroup[i] = 0;
        BigThanLB[i] = 0;
        groupElement[i] = i;
        tmpUB[i] = UB[i];
        scSizeGroup[i] = 0;
    }

    // Initialize s1 with partition1
    for (i = 0; i < N; i++) {
        s1[i] = partition1[i];
    }
    BuildDeltaMatrix(s1);
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            Delta_Matrix_p1[i][j] = Delta_Matrix[i][j];
        }
    }
    BuildGroupDiversityForCrossover(s1, groupDiversity_s1);

    // Initialize s2 with partition2
    for (i = 0; i < N; i++) {
        s2[i] = partition2[i];
    }
    BuildDeltaMatrix(s2);
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            Delta_Matrix_p2[i][j] = Delta_Matrix[i][j];
        }
    }
    BuildGroupDiversityForCrossover(s2, groupDiversity_s2);

    int targetGroup = -1;
    // Main crossover process
    for (i = 0; i < K; i++) {
        elementCount = 0;
        if (uniform_rnd_number() < 0.5) {
            // Process partition 1
            process_partition(groupDiversity_s1, s1, tmpUB, childSolution, vectorElement, K, N, &elementCount, &targetGroup);  
        } else {
            // Process partition2 (similar to partition1 logic)
            process_partition(groupDiversity_s2, s2, tmpUB, childSolution, vectorElement, K, N, &elementCount, &targetGroup);  
        }

        // Update group diversity
        for (j = 0; j < elementCount; j++) {
            groupDiversity_s1[s1[tmpEle[j]]] -= Delta_Matrix_p1[tmpEle[j]][s1[tmpEle[j]]];
            groupDiversity_s2[s2[tmpEle[j]]] -= Delta_Matrix_p2[tmpEle[j]][s2[tmpEle[j]]];
            s1[tmpEle[j]] = -1;
            s2[tmpEle[j]] = -1;
        }

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

    // Assign unprocessed elements to meet lower bounds Pseudo code line 22
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
            if (childSolution[j] == targetGroup) {
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = random_int(elementCount);
        childSolution[SelectEle[selectedElement]] = -1;
        vectorElement[SelectEle[selectedElement]] = SelectEle[selectedElement];
        scSizeGroup[targetGroup]--;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            BigThanLB[targetGroup] = 0;
        }
        processedCount++;
    }

    // Assign elements to meet lower bounds Pseudo code line 22
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
        childSolution[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            LBGroup[targetGroup] = 0;
        }
        totalBelowLowerBound++;
    }

    // Assign elements to meet upper bounds Pseudo code line 23
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
        childSolution[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == UB[targetGroup]) {
            UBGroup[targetGroup] = 0;
        }
        totalSize++;
    }
}

double LocalSearchCriterionCalculation(Solution* sol1, Solution* sol2) {
    /*
     * Evaluates the quality and dissimilarity of partitions.
     * It calculates the value that combines the ratio of costs (sol1->f / sol2->f)
     * and a dissimilarity factor between sol1->partition and sol2->partition.
     */

    // Handle potential division by zero
    if (sol2->objective == 0.0) {
        Rprintf("Error: Division by zero (sol2->f is zero).\n");
        return -1;
    }

    int i, j;
    int totalPairs = (N * (N - 1)) / 2;  // Number of unique pairs (i, j) where i < j
    int count = 0;

    // Loop over all pairs of elements to count dissimilar pairs
    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            // Count cases where elements are grouped differently in the two partitions
            if ((sol1->s[i] == sol1->s[j]) != (sol2->s[i] == sol2->s[j])) {
                count++;
            }
        }
    }

    // Calculate dissimilarity factor
    double dissimilarityFactor = ((double)count / totalPairs) * K;

    // Calculate and return the criterion value (ratio of costs + weighted dissimilarity factor)
    return sol1->objective / sol2->objective + alpha * dissimilarityFactor;
}

void ClearDeltaMatrix(void) {
	/* Resets the delta_f matrix */
    for (int i = 0; i < N; ++i) {
        // should this not be over N ?!
        for (int j = 0; j < K; j++) {
            Delta_Matrix[i][j] = 0.0;
        }
    }
}

void BuildDeltaMatrix(int partition[]) {
	/*  Builds the delta_f matrix and calculates the objective function value */

	ClearDeltaMatrix();

    int i, j;
	// Update Delta_Matrix based on distances
    for (int i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Delta_Matrix[i][partition[j]] += Distances[i][j];
        }
    }

    // Calculate the objective function value
    f_objective = 0.0;
    for (i = 0; i < N; i++) {
        f_objective += Delta_Matrix[i][partition[i]];
    }
    f_objective /= 2.0;
}

void BuildGroupDiversityForCrossover(int partition[], double groupDiversity[]) {
	/*  Builds group diversity values for crossover */
	
    // Initialize group diversity values to zero
	for (int i = 0; i < K; i++) groupDiversity[i] = 0.0;
	
	// Compute group diversity based on distances
    for (int i = 0; i < N; i++) {
        int group_i = partition[i];
        for (int j = 0; j < N; j++) {
            if (group_i == partition[j]) {
                groupDiversity[group_i] += Distances[i][j];
            }
        }
    }
}

void OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup) {
    // Update Delta_Matrix for all elements affected by the move
	for (int j = 0; j < N; j++) {
		if (j != i) {
			Delta_Matrix[j][oldGroup] -= Distances[i][j];
			Delta_Matrix[j][newGroup] += Distances[i][j];
		}
	}
}

void AssignMemoryDiversity(void) {
    /*  Allocates memory dynamically for various arrays and matrices necessary 
	for the algorithm's execution. This includes structures for population management, 
	distance matrices, diversity measures, and neighborhood exploration.
	*/
    
    S = (Solution*)malloc(beta_max * sizeof(Solution));
    O = (Solution*)malloc(beta_max * sizeof(Solution));
    int i;
    for (i = 0; i < beta_max; i++) {
        S[i].s = (int*)malloc(N * sizeof(int));
        O[i].s = (int*)malloc(N * sizeof(int));
        S[i].SizeG = (int*)malloc(K * sizeof(int));
        O[i].SizeG = (int*)malloc(K * sizeof(int));
    }    
    
    Delta_Matrix = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix_p1 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p1[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix_p2 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p2[i] = (double*)malloc(K * sizeof(double));
    groupDiversity_s1 = (double*)malloc(K * sizeof(double));
    groupDiversity_s2 = (double*)malloc(K * sizeof(double));
    

    S_best.s = (int*)malloc(N * sizeof(int));
    S_best.SizeG = (int*)malloc(K * sizeof(int));
    
    Avg = (double**)malloc(K * sizeof(double*));
    for (i = 0; i < K; i++) Avg[i] = (double*)malloc(K * sizeof(double));
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

void ReleaseMemoryDiversity(void) {
    /* responsible for reading the input file, 
    initializing matrices, and setting constraints on group sizes. */ 
    
    int i;
    for (i = 0; i < beta_max; i++) {
        free(S[i].s); S[i].s = NULL;
        free(S[i].SizeG); S[i].SizeG = NULL;
        free(O[i].s); O[i].s = NULL;
        free(O[i].SizeG); O[i].SizeG = NULL;
    }
    free(S); S = NULL;
    free(O); O = NULL;
    
    free(S_best.s); S_best.s = NULL;
    free(S_best.SizeG); S_best.SizeG = NULL;
    free(LB); LB = NULL;
    free(UB); UB = NULL;
    
    for (i = 0; i < N; i++) {
        free(Delta_Matrix[i]); Delta_Matrix[i] = NULL;
        free(Delta_Matrix_p1[i]); Delta_Matrix_p1[i] = NULL;
        free(Delta_Matrix_p2[i]); Delta_Matrix_p2[i] = NULL;
    }
    free(Delta_Matrix); Delta_Matrix = NULL;
    free(Delta_Matrix_p1); Delta_Matrix_p1 = NULL;
    free(Delta_Matrix_p2); Delta_Matrix_p2 = NULL;
    free(groupDiversity_s1); groupDiversity_s1 = NULL;
    free(groupDiversity_s2); groupDiversity_s2 = NULL;
    free(Avg); Avg = NULL;
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

// Generate a random integer from zero to max-1
int random_int(int max) {
  GetRNGstate();
  double my_number = unif_rand();
  PutRNGstate();
  return (int) floor(my_number * max);
}
