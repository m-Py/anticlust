#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "fifr_header.h"
#include <R.h>
#include <Rinternals.h>

// Global variables (from the paper)
static int beta_max;
static int N;
int M; // elements, groups
static double** Distances;
static double** DistancesT;
static int* LB; // lower bound
static int* UB; // upper bound
static double theta, theta_max, theta_min;
static double alpha;
static int beta_min;
double phi;
int tau;
int L;

// Matrices
static double **Delta_Matrix;
static double **Delta_Matrix_p1;
static double **Delta_Matrix_p2;

// Solutions
static Solution S_best;
static Solution *S;
static Solution *O;

// Other variables
static double f_objective;
int noImpCounter;
int kmax;
static int maxNumberIterations;
static double Time_limit;

// For crossover
static int *vectorElement;
static int *groupElement;
static int *SelectGroup;
static int *SelectEle;
static int *tmpEle;
static int *s1;
static int *s2;
static double *groupDiversity_s1;
static double *groupDiversity_s2;
static int *LBGroup, *UBGroup, *tmpUB, *BigThanLB;

// IFR search
int knn;
int *p;

void FIFR_ClearDeltaMatrix(void);
void FIFR_BuildDeltaMatrix(int partition[]);
void FIFR_OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup);
void FIFR_BuildGroupDiversityForCrossover(int partition[], double groupDiversity[]);
void RandlsDiversity(int s[], int SizeGroup[], double *objective);
void FIFR_CrossoverDiversity(int partition1[], int partition2[], int childSolution[], int scSizeGroup[]);
void FIFR_process_partition(double* groupDiversity, int* partition, int* tmpUB, int* childSolution, int* vectorElement,
    int M, int N, int *element_count, int *target_group);
void BreakGroupConstraints(int partition[], int SizeGroup[], double *objective);
void FitGroupConstraints(int partition[], int SizeGroup[], double *objective);
void FIFR_SearchAlgorithmDiversity(void);
void FIFR_AssignMemoryDiversity(void);
void FIFR_ReleaseMemoryDiversity(void);

/* FIFR for Anticlustering based on a Distance matrix
 * 
 * param *distannces: vector of data points (in R, this is a distance matrix,
 *         the matrix structure must be restored in C)
 * param *N_in: The number of elements (i.e., number of "rows" in *data)
 * param *M_in: The number of clusters. When lower_bound and upper_boound are set to the number of M, 
 *              the clusters will be equally sized. 
 * param *number_of_iterations: A number that defines how many times the steps in the search algorithm are repeated.
 * param *clusters: A predefined vector of length M specifies the number of elements in each cluster.
 *               If a default vector [-1] is provided, cluster sizes will be determined based 
 *               on the lower and upper bounds. When a cluster size array is provided, 
 *               the lower and upper bounds are ignored as they become redundant.
 * param *lower_bound: Minimum number of elements in each anticluster. 
 * param *upper_bound: Maximum number of elements in each anticluster.
 * param *Beta_max: The algorithm begins with a pool of random initial solutions of size beta_max. 
 *                   Over time, the size of the solution pool decreases linearly until it reaches beta_min.
 * param *Theta_max: Parameter for the strength of undirected perturbation,
 *                   which decreases linearly over time from theta_max to theta_min..
 * param *Theta_min: Parameter for the strength of undirected perturbation, 
 *                   which decreases linearly over time from theta_max to theta_min..
 * param *Beta_min: The minimum solution pool size the algorithm should reach before making a determination.
 * param *Phi: A parameter that determines the population size when initiating a new round of exploration upon
 *              triggering the jump-back mechanism
 * param *Tau: A parameter when "noImp" exceeds it, the algorithm transits from exploitation to exploration strategy
 * param *Kmax: A parameter that determines the maximum degree of constraint violation in the IFR search
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

void feasible_and_infeasible_region_search(
    double *distances,
    int *N_in,
    int *M_in,
    int *number_of_iterations,
    int *clusters,
    int *upper_bound, 
    int *lower_bound, 
	int *Beta_max, 
	int *elapsed_time,
	double *Theta_max,
	double *Theta_min,
	int *Beta_min,
	double *Phi,
    int *Tau,
    int *Kmax,
    double *Alpha,
	int *result,
	double *score,
	int *mem_error){

    N = *N_in;
    M = *M_in;
    beta_max = *Beta_max;
    theta = *Theta_max;
    theta_max = *Theta_max;
    theta_min = *Theta_min;
    beta_min = *Beta_min;
    phi = *Phi;
    tau = *Tau;
    Time_limit = *elapsed_time;
    kmax = *Kmax;
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
        LB = (int*)malloc(M * sizeof(int));
        if (LB == NULL) { *mem_error = 1; return; }
        UB = (int*)malloc(M * sizeof(int));
        if (UB == NULL) { *mem_error = 1; return; }
        for (int i = 0; i < M; i++) {
            LB[i] = *lower_bound;  // Assuming lower_bound is a pointer to an int
            UB[i] = *upper_bound;  // Assuming upper_bound is a pointer to an int
        }
    } else {
        LB = (int*)malloc(M * sizeof(int));
        if (LB == NULL) { *mem_error = 1; return; }
        UB = (int*)malloc(M * sizeof(int));
        if (UB == NULL) { *mem_error = 1; return; }
        for (int i = 0; i < M; i++) {
            LB[i] = clusters[i];  // Assuming lower_bound is a pointer to an int
            UB[i] = clusters[i];  // Assuming upper_bound is a pointer to an int
        }
    }
    
    FIFR_AssignMemoryDiversity();
    if (*mem_error == 1) {
        return;
    }
        
    FIFR_SearchAlgorithmDiversity();
    
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
    
    FIFR_ReleaseMemoryDiversity();
}

void FIFR_SearchAlgorithmDiversity() {

    /* Algorithm 1: Main procedure of FIFR. */

    int beta = beta_max;
    int i,j,k;
    int pickedSolution;

    //important! for windows and linux there is a differnt definition of this time
    //on windows its the wall time, on linux the CPU time
    S_best.objective = -INFINITY;
    
    // Initial population generation
    for (int i = 0; i < beta; i++) {
        /* Algorithm 2: random initial solution construction*/
        InitialSolution(S[i].s, S[i].SizeG);
        RandlsDiversity(S[i].s, S[i].SizeG, &(S[i].objective));
        
        if (S_best.objective < S[i].objective) {
            for (int j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
            for (int k = 0; k < M; k++) S_best.SizeG[k] = S[i].SizeG[k];
            S_best.objective = S[i].objective;
        }
        
    }
    
    noImpCounter = 0;
    for (int counter = 1; counter <= maxNumberIterations; counter++) {
        L = (int)(theta * N / M);

        // Exploration
        if (beta > beta_min) {
            /* Algorithm 4-1: Perturbation*/
            for (int i = 0; i < beta; i++) {
                Perturbation(L, S[i].s, S[i].SizeG);
                RandlsDiversity(S[i].s, S[i].SizeG, &(S[i].objective));
                
                if (S_best.objective < S[i].objective) {
                    for (int j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                    for (int k = 0; k < M; k++) S_best.SizeG[k] = S[i].SizeG[k];
                    S_best.objective = S[i].objective;
                }
            }
            
            /* Algorithm 4-2: Crossover*/
            if (beta > 1){
                for (i = 0; i < beta; i++){
                    pickedSolution = FIFR_random_int(beta);
                    if (pickedSolution == i){
                        pickedSolution = (pickedSolution + 1) % beta;
                    }

                    FIFR_CrossoverDiversity(S[i].s, S[pickedSolution].s, O[i].s, O[i].SizeG);
                    RandlsDiversity(O[i].s, O[i].SizeG, &(O[i].objective));
                }
                for (i = 0; i < beta; i++){
                    if (O[i].objective >= S[i].objective){
                        for (j = 0; j < N; j++) S[i].s[j] = O[i].s[j];
                        for (k = 0; k < M; k++) S[i].SizeG[k] = O[i].SizeG[k];
                        S[i].objective = O[i].objective;
                    } else if (FIFR_LocalSearchCriterionCalculation(&O[i], &S[i]) > 1) {
                        for (j = 0; j < N; j++) S[i].s[j] = O[i].s[j];
                        for (k = 0; k < M; k++) S[i].SizeG[k] = O[i].SizeG[k];
                        S[i].objective = O[i].objective;
                    }

                    if (S_best.objective < S[i].objective){
                        for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                        for (k = 0; k < M; k++) S_best.SizeG[k] = S[i].SizeG[k];
                        S_best.objective = S[i].objective;
                    }
                }
            }
        }

        /* Algorithm 5: Exploitation */
        while (beta == beta_min){
            for (i = 0; i < beta; i++){
                Perturbation(L, S[i].s, S[i].SizeG);
                RandlsDiversity(S[i].s, S[i].SizeG, &(S[i].objective));

                if (S_best.objective < S[i].objective){
                    for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                    for (k = 0; k < M; k++) S_best.SizeG[k] = S[i].SizeG[k];
                    S_best.objective = S[i].objective;
                    noImpCounter = 0;
                } else if (S_best.objective >= S[i].objective) {
                    noImpCounter++;
                }
            }

            if (noImpCounter > tau) {
                beta = (int) (phi * beta_max);
                noImpCounter = 0;
                break;
            }
        }
        
        /* Algorithm 6: Infeasible Region Search */
        for (i = 0; i < beta; i++) {
            knn = 1;
            if(knn < kmax){
                BreakGroupConstraints(S[i].s, S[i].SizeG, &(S[i].objective));
                FitGroupConstraints(S[i].s, S[i].SizeG, &(S[i].objective));
                if (S_best.objective < S[i].objective) {
                    for (j = 0; j < N; j++) S_best.s[j] = S[i].s[j];
                    for (k = 0; k < M; k++) S_best.SizeG[k] = S[i].SizeG[k];
                    S_best.objective = S[i].objective;
                    knn = kmax;
                } else {
                    knn++;
                }
            }
        }
        qsort(S, beta, sizeof(Solution), FIFR_CompareSolution);
        //int previousBeta = beta;
        beta = (int)(beta - (beta - 1) * counter / maxNumberIterations);
        theta = theta_max - (theta_max - theta_min) * counter / maxNumberIterations;
    }
}

void FIFR_swap_elements(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void FIFR_fisher_yates_shuffle(int arr[], int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = FIFR_random_int(i + 1);
        FIFR_swap_elements(&arr[i], &arr[j]);
    }
}

void InitialSolution(int s[], int SizeG[]){
    /* Algorithm 2: Initial random solution generation that respects the group size constraints */

    // Allocate memory
    int *isAssigned = (int *)malloc(N * sizeof(int)); // Tracks if an element is assigned
    int *groupSize = (int *)malloc(M * sizeof(int)); // Stores the size of each group
    int *permutedIndexList = (int *)malloc(N * sizeof(int)); // Stores the size of each group
    int *permutedGroupList = (int *)malloc(M * sizeof(int)); // Stores the size of each group

    int i;
    for (i = 0; i < M; i++) groupSize[i] = 0;
    for (i = 0; i < N; i++) isAssigned[i] = 0;
    for (i = 0; i < N; i++) permutedIndexList[i] = i;
    for (i = 0; i < M; i++) permutedGroupList[i] = i;

    FIFR_fisher_yates_shuffle(permutedIndexList, N);
    FIFR_fisher_yates_shuffle(permutedGroupList, M);

    // calculate the total number of elements that need to satisfy the lower bounds
    int total_assigned = 0;
    int total_LB = 0;
    for (i = 0; i < M; i++) {
        total_LB += LB[i];
    }

    // First phase: Assign elements to satisfy lower bound constraints (LB)
    int selected_element = 0;
    while (total_assigned < total_LB){
        for (int group = 0; group < M; group++){
            if (groupSize[group] < LB[group]){
                s[permutedIndexList[selected_element]] = group;
                isAssigned[permutedIndexList[selected_element]] = 1;
                groupSize[group]++;
                total_assigned++;
                break; // Move to the next element once assigned
            }
        }
        selected_element++;
    }

    // Second phase: Assign remaining elements randomly while respecting upper bounds (UB)
    while (total_assigned < N){
        for (int group = 0; group < M; group++) {
            if (groupSize[permutedGroupList[group]] < UB[permutedGroupList[group]]) {
                s[permutedIndexList[selected_element]] = permutedGroupList[group];
                isAssigned[permutedIndexList[selected_element]] = 1;
                groupSize[permutedGroupList[group]]++;
                total_assigned++;
                break;
            }
        }
        FIFR_fisher_yates_shuffle(permutedGroupList, M); // Shuffle groups for randomness
        selected_element++;
    }

    // Copy final group sizes into the output array SizeG
    for (i = 0; i < M; i++) SizeG[i] = groupSize[i];

    // Free allocated memory
    free(groupSize); groupSize = NULL;
    free(isAssigned); isAssigned = NULL;
    free(permutedIndexList); permutedIndexList = NULL;
    free(permutedGroupList); permutedGroupList = NULL;
}

void RandlsDiversity(int s[], int SizeGroup[], double *objective){
    /* Algorithm 3: Neighborhood search to find local optima */
    const double DELTA_THRESHOLD = 0.0001; // Define a constant for comparison threshold
    int v, g, u;
    int oldGroup, oldGroup1, t;

    // Build the delta_f matrix for objetive changes
    FIFR_BuildDeltaMatrix(s);

    // Initialize the delta_f value
    double delta_f = -99999.0;

    int imp;
    do {
        imp = 0; // Reset improvement flag

        // First loop: Move individual elements to improve partition
        for (v = 0; v < N; v++){
            for (g = 0; g < M; g++){
                if ((s[v] != g) && (SizeGroup[s[v]] > LB[s[v]]) && (SizeGroup[g] < UB[g])) {
                    delta_f = Delta_Matrix[v][g] - Delta_Matrix[v][s[v]];
                    if (delta_f > DELTA_THRESHOLD) {
                        oldGroup = s[v];

                        // Update delta_f matrix for the move
                        FIFR_OneMoveUpdateDeltaMatrix(v, oldGroup, g);

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
        for (v = 0; v < N; v++){
            for(u = v + 1; u < N; u++){
                // Only swap if nodes are in differnt groups
                if (s[v] != s[u]){
                    delta_f = (Delta_Matrix[v][s[u]] - Delta_Matrix[v][s[v]]) 
                    + (Delta_Matrix[u][s[v]] - Delta_Matrix[u][s[u]]) 
                    - DistancesT[v][u];

                    if (delta_f > DELTA_THRESHOLD){
                        oldGroup = s[v];
                        oldGroup1 = s[u];

                        // Update delta_f matrix for the swap
                        FIFR_OneMoveUpdateDeltaMatrix(v, oldGroup, oldGroup1);
                        FIFR_OneMoveUpdateDeltaMatrix(u, oldGroup1, oldGroup);

                        // Swap the elements in the partition
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
    } while (imp == 1); // Continue until no improvement is made

    // Update partition array with the final assignments
    *objective = f_objective;
}

void Perturbation(int L, int s[], int SizeGroup[]){

    int perturb_type;
    int v, g, x, y;
    int NumberNeighbors, oldGroup, swap;
    int theta, count = 0;

    theta = L;
    NumberNeighbors = N * (N - 1) / 2 + N * M;

    while (count < theta){
        perturb_type = FIFR_random_int(NumberNeighbors);

        if (perturb_type < N * M){ // Type 1: Random (element, group) perturbation
            v = FIFR_random_int(N); // Randomly choose an element v
            g = FIFR_random_int(M); // Randomly choose a group g

            if (s[v] != g && SizeGroup[s[v]] > LB[s[v]] && SizeGroup[g] < UB[g]){
                oldGroup = s[v];
                SizeGroup[oldGroup] -= 1;
                SizeGroup[g] += 1;
                s[v] = g;
                count++;
            }
        }
        else { // Type 2: Random (element, element) swap perturbation
            x = FIFR_random_int(N); // Randomly choose element x
            y = FIFR_random_int(N); // Randomly choose element y

            // Apply partition if elements are in different groups
            if (s[x] != s[y] && x != y){
                swap = s[x];
                s[x] = s[y];
                s[y] = swap;
                count++;
            }
        }
    }
}

void FIFR_CrossoverDiversity(int partition1[], int partition2[], int childSolution[], int scSizeGroup[]){
    int i, j;
    int elementCount, processedCount, selectedElement;
    int totalLowerBound, totalBelowLowerBound;

    // Initialize arrays
    for (i = 0; i < N; i++){
        vectorElement[i] = i;
        childSolution[i] = -1;
    }
    for (i = 0; i < M; i++){
        LBGroup[i] = 0;
        UBGroup[i] = 0;
        BigThanLB[i] = 0;
        groupElement[i] = 0;
        tmpUB[i] = UB[i];
        scSizeGroup[i] = 0;
    }

    // Intialize s1 with partition1
    for (i = 0; i < N; i++) {
        s1[i] = partition1[i];
    }
    FIFR_BuildDeltaMatrix(s1);
    for (i = 0; i < N; i++){
        for (j = 0; j < M; j++){
            Delta_Matrix_p1[i][j] = Delta_Matrix[i][j];
        }
    }

    FIFR_BuildGroupDiversityForCrossover(s1, groupDiversity_s1);

    // Intialize s2 with partition2
    for (i = 0; i < N; i++) {
        s2[i] = partition2[i];
    }
    FIFR_BuildDeltaMatrix(s2);
    for (i = 0; i < N; i++){
        for (j = 0; j < M; j++){
            Delta_Matrix_p2[i][j] = Delta_Matrix[i][j];
        }
    }

    FIFR_BuildGroupDiversityForCrossover(s2, groupDiversity_s2);

    int targetGroup = -1;
    // Main crossover process
    for (i = 0; i < M; i++){
        elementCount = 0;
        if (uniform_rnd_number() < 0.5){
            // process partition1
            FIFR_process_partition(groupDiversity_s1, s1, tmpUB, childSolution, vectorElement, M, N, &elementCount, &targetGroup);
        } else {
            // process partition2
            FIFR_process_partition(groupDiversity_s2, s2, tmpUB, childSolution, vectorElement, M, N, &elementCount, &targetGroup);
        }

        // Update group diversity
        for (j = 0; j < elementCount; j++){
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
    for (i = 0; i < M; i++){
        totalLowerBound += LB[i];
        if (scSizeGroup[i] < LB[i]){
            processedCount += scSizeGroup[i];
            totalBelowLowerBound += scSizeGroup[i];
            LBGroup[i] = 1;
        } else {
            processedCount += LB[i];
        }
        if (scSizeGroup[i] > LB[i]){
            BigThanLB[i] = 1;
        }
    }

    // Assign unprocessed elements to meet lower bounds
    for (i = 0; i < N; i++){
        if (vectorElement[i] != -1){
            processedCount++;
        }
    }
    while (processedCount < totalLowerBound){
        targetGroup = FIFR_random_int(M);
        while (BigThanLB[targetGroup] == 0){
            targetGroup = (targetGroup + 1) % M;
        }

        elementCount = 0;
        for (j = 0; j < N; j++){
            if (childSolution[j] == targetGroup){
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = FIFR_random_int(elementCount);
        childSolution[SelectEle[selectedElement]] = -1;
        vectorElement[SelectEle[selectedElement]] = SelectEle[selectedElement];
        scSizeGroup[targetGroup]--;
        if (scSizeGroup[targetGroup] == LB[targetGroup]){
            BigThanLB[targetGroup] = 0;
        }
        processedCount++;
    }

    // Assign elements to meet lower bounds
    int sumLB = 0;
    for (i = 0; i < M; i++){
        if (LBGroup[i] == 1){
            sumLB += LB[i];
        }
    }
    while (totalBelowLowerBound < sumLB){
        targetGroup = FIFR_random_int(M);
        while (LBGroup[targetGroup] == 0){
            targetGroup = (targetGroup + 1) % M;
        }

        elementCount = 0;
        for (i = 0; i < N; i++){
            if (vectorElement[i] != -1){
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = FIFR_random_int(elementCount);
        childSolution[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == LB[targetGroup]){
            LBGroup[targetGroup] = 0;
        }
        totalBelowLowerBound++;
    }

    // Assign elements to meet upper bounds
    int totalSize = 0;
    for (i = 0; i < M; i++){
        totalSize += scSizeGroup[i];
        if (scSizeGroup[i] < UB[i]){
            UBGroup[i] = 1;
        }
    }
    while (totalSize < N){
        targetGroup = FIFR_random_int(M);
        while (UBGroup[targetGroup] == 0){
            targetGroup = (targetGroup + 1) % M;
        }

        elementCount = 0;
        for (i = 0; i < N; i++){
            if (vectorElement[i] != -1){
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = FIFR_random_int(elementCount);
        childSolution[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == UB[targetGroup]){
            UBGroup[targetGroup] = 0;
        }
        totalSize++;
    }
}

void FIFR_BuildGroupDiversityForCrossover(int partition[], double groupDiversity[]){
    /* Builds the group diversity values for crossover */
    int i, j, group_i;
    // Initialize group diversity values to zero
    for (i = 0; i < M; i++) groupDiversity[i] = 0.0;

    for (i=0; i<N; i++){
        group_i = partition[i];
        for (j=i+1; j<N; j++){
            if (group_i == partition[j]) {
                groupDiversity[partition[i]] += Distances[i][j];
            }
        }
    }
}

void FIFR_process_partition(
    double* groupDiversity,
    int* partition,
    int* tmpUB,
    int* childSolution,
    int* vectorElement,
    int M,
    int N,
    int *element_count,
    int *target_group
) {
    int i, selectedGroup, processedCount, selectedElement;

    int elementCount = *element_count;
    int targetGroup = *target_group;
    double maxGroupDiversity = -1e100;
    for (i = 0; i < M; i++){
        if (groupDiversity[i] > maxGroupDiversity){
            maxGroupDiversity = groupDiversity[i];
            selectedGroup = i;
        }
    }

    for (i = 0; i < N; i++){
        if (partition[i] == selectedGroup){
            SelectEle[elementCount++] = i;
        }
    }

    int groupCount = 0;
    for(i = 0; i < M; i++){
        if(tmpUB[i] != -1 && tmpUB[i] >= elementCount){
            SelectGroup[groupCount++] = i;
        }
    }

    if (groupCount == 0){ // no valid group found
        int minDiff = 999999;
        for (i = 0; i < M; i++){
            if (tmpUB[i] != -1 && elementCount - tmpUB[i] < minDiff){
                minDiff = elementCount - tmpUB[i];
                targetGroup = i;
            }
        }

        processedCount = 0;
        while (processedCount < elementCount - minDiff){
            selectedElement = FIFR_random_int(elementCount);
            while (SelectEle[selectedElement] == -1){
                selectedElement = (selectedElement + 1) % elementCount;
            }

            childSolution[SelectEle[selectedElement]] = targetGroup;
            tmpEle[processedCount++] = SelectEle[selectedElement];
            vectorElement[SelectEle[selectedElement]] = -1;
            SelectEle[selectedElement] = -1;
        }
        elementCount = processedCount;
    } else {
        targetGroup = SelectGroup[FIFR_random_int(groupCount)];
        for (i = 0; i < elementCount; i++){
            childSolution[SelectEle[i]] = targetGroup;
            vectorElement[SelectEle[i]] = -1;
            tmpEle[i] = SelectEle[i];
        }
    }
    *element_count = elementCount;
    *target_group = targetGroup;
}

double FIFR_LocalSearchCriterionCalculation(Solution* sol1, Solution* sol2){
    /*
     * Evaluates the quality and dissimilarity of partitions.
     * It calculates the value that combines the ratio of costs (sol1->f / sol2->f)
     * and a dissimilarity factor between sol1->partition and sol2->partition.
    */

    // Handle potential division by zero
    if (sol2->objective == 0) {
        return -1;
    }

    int i, j;
    int totalPairs = N * (N - 1) / 2; // Number of unique paris (i, j) with i < j
    int count = 0;

    // Loop over all pairs of elements to count dissimilar pairs
    for (i = 0; i < N - 1; i++){
        for (j = i + 1; j < N; j++){
            if ((sol1->s[i] == sol1->s[j]) != (sol2->s[i] == sol2->s[j])){
                count++;
            }
        }
    }

    // Calculate dissimilarity factor
    double dissimilarityFactor = ((double)count / totalPairs) * M;

    // Calculate and return the final criterion value (ratio if costs + weighted dissimilarity factor)
    return sol1->objective / sol2->objective + alpha * dissimilarityFactor;
}

void BreakGroupConstraints(int partition[], int SizeGroup[], double *objective){
    int i, v, g;
    int oldGroup;
    double delta_f = -999999.0;
    int imp;
    for (i = 0; i < N; i++) p[i] = partition[i];
    FIFR_BuildDeltaMatrix(p);

    do {
        imp = 0;
        for(v = 0; v < N; v++){
            for(g = 0; g < M; g++){
                if( (p[v] != g) && (SizeGroup[p[v]] > (LB[p[v]] - knn)) && (SizeGroup[g] < (UB[g] + knn))) {
                    delta_f = Delta_Matrix[v][g] - Delta_Matrix[v][p[v]];
                    if(delta_f > 0.0001){
                        oldGroup = p[v];
                        FIFR_OneMoveUpdateDeltaMatrix(v, oldGroup, g);
                        SizeGroup[oldGroup] = SizeGroup[oldGroup] - 1;
                        SizeGroup[g] = SizeGroup[g] + 1;
                        p[v] = g;
                        f_objective += delta_f;
                        imp = 1;
                    }
                }
            }    
        }
    } while (imp == 1);
    FIFR_BuildDeltaMatrix(p);
    *objective = f_objective;
    for (i = 0; i < N; i++) partition[i] = p[i];
}

void FitGroupConstraints(int partition[], int SizeGroup[], double *objective){
    int i,j;

    for (i = 0; i < N; i++) p[i] = partition[i];
    FIFR_BuildDeltaMatrix(p);

    // Violation-Arrays
    int *needLB = (int*)calloc(M, sizeof(int));
    int *excessUB = (int*)calloc(M, sizeof(int));
    if (!needLB || !excessUB) { perror("calloc"); return; }

    int numa = 0; // missing LB-places
    int numb = 0; // missing UB-excesses

    for (i = 0; i < M; i++) {
        if (SizeGroup[i] < LB[i]) {
            needLB[i] = LB[i] - SizeGroup[i];
            numa += needLB[i];
        } else if (SizeGroup[i] > UB[i]) {
            excessUB[i] = SizeGroup[i] - UB[i];
            numb += excessUB[i];
        }
    }

    // collect element from groups with UB + knn > UB
    int up_elements = 0;
    for (i = 0; i < M; i++) 
        if (excessUB[i] > 0) up_elements += SizeGroup[i];

    int *gene = NULL;
    if (up_elements > 0){
        gene = (int*)malloc(up_elements * sizeof(int));
        if (!gene) { perror("malloc"); return; }
        int index = 0;
        for (i = 0; i < M; i++){
            if (excessUB[i] > 0){
                for (j = 0; j < N; j++){
                    if (p[j] == i) gene[index++] = j;
                }
            }
        }
    }

    // Step 1: 
    while (numb > 0){
        int kn = -1; // needing element
        int kg = -1; // target group
        double bestMin = 1e100;

        // find element which least contribution to group diversity
        for (i = 0; i < up_elements; i++){
            int v = gene[i];
            if (v == -1) continue;
            int gv = p[v];
            if (excessUB[gv] <= 0) continue;

            // Schutz gegen div0 (sollte praktisch nicht passieren)
            double denom = (SizeGroup[gv] > 0) ? (double)SizeGroup[gv] : 1.0;
            double score = Delta_Matrix[v][gv] / denom;

            if (score < bestMin) {
                bestMin = score;
                kn = v;
            }
        }

        if (kn == -1) break;

        // kn aus gene "entfernen"
        for (i = 0; i < up_elements; i++) if (gene[i] == kn) { gene[i] = -1; break; }

        // best target group 
        double bestMax = -1e100;
        for (j = 0; j < M; j++) {
            if (SizeGroup[j] >= UB[j]) continue;

            double denom = (SizeGroup[j] > 0) ? (double)SizeGroup[j] : 1.0;
            double score = Delta_Matrix[kn][j] / denom;

            if (score > bestMax) {
                bestMax = score;
                kg = j;
            }
        }

        if (kg == -1) break;

        int oldg = p[kn];

        // move kn: oldg -> kg
        FIFR_OneMoveUpdateDeltaMatrix(kn, oldg, kg);
        p[kn] = kg;

        SizeGroup[oldg]--;
        SizeGroup[kg]++;

        excessUB[oldg]--; 
        numb--;

        // if target group was under LB, this also reduces LB-deficit
        if (needLB[kg] > 0) { needLB[kg]--; numa--; }
    }

    // Step 2: LB auffüllen
    while (numa > 0)
    {
        int underg = -1;
        for (i = 0; i < M; i++) {
            if (needLB[i] > 0) { underg = i; break; }
        }
        if (underg == -1) break;

        int chosei = -1;
        double bestDelta = -1e100;

        // search for an element that can move to underg and has the best improvement
        for (j = 0; j < N; j++) {
            int gj = p[j];
            if (gj == underg) continue;

            if (SizeGroup[gj] > LB[gj] && SizeGroup[gj] <= UB[gj]) {
                double delt = Delta_Matrix[j][underg] - Delta_Matrix[j][gj];
                if (delt > bestDelta) {
                    bestDelta = delt;
                    chosei = j;
                }
            }
        }

        if (chosei == -1) break;

        int oldg = p[chosei];

        FIFR_OneMoveUpdateDeltaMatrix(chosei, oldg, underg);
        p[chosei] = underg;

        SizeGroup[oldg]--;
        SizeGroup[underg]++;

        needLB[underg]--;
        numa--;
    }

    // calculate final objective
    FIFR_BuildDeltaMatrix(p);
    *objective = f_objective;

    // copy back
    for (i = 0; i < N; i++) partition[i] = p[i];

    free(needLB);
    free(excessUB);
    free(gene);
}

/* 
    -------------------------------------------
    ------------- Helper function -------------
    ------------------------------------------- 
*/

int FIFR_CompareSolution(const void *first, const void *second) {
    Solution *solution1 = (Solution*)first;
    Solution *solution2 = (Solution*)second;
    
    if (solution1->objective < solution2->objective) return 1;
    else if (solution1->objective > solution2->objective) return -1;
    else return 0;
}

void FIFR_OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup) {
    for (int j = 0; j < N; j++) {
        if (j != i) {
            Delta_Matrix[j][oldGroup] -= Distances[i][j];
            Delta_Matrix[j][newGroup] += Distances[i][j];
        }
    }
}

void FIFR_ClearDeltaMatrix() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            Delta_Matrix[i][j] = 0.0;
        }
    }
}

void FIFR_BuildDeltaMatrix(int partition[]) {
    FIFR_ClearDeltaMatrix();
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Delta_Matrix[i][partition[j]] += Distances[i][j];
        }
    }
    
    f_objective = 0.0;
    for (int i = 0; i < N; i++) {
        f_objective += Delta_Matrix[i][partition[i]];
    }
    f_objective = f_objective / 2.0;
}

void FIFR_AssignMemoryDiversity(void) {
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
        S[i].SizeG = (int*)malloc(M * sizeof(int));
        O[i].SizeG = (int*)malloc(M * sizeof(int));
    }    
    
    Delta_Matrix = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix[i] = (double*)malloc(M * sizeof(double));
    Delta_Matrix_p1 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p1[i] = (double*)malloc(M * sizeof(double));
    Delta_Matrix_p2 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p2[i] = (double*)malloc(M * sizeof(double));
    groupDiversity_s1 = (double*)malloc(M * sizeof(double));
    groupDiversity_s2 = (double*)malloc(M * sizeof(double));
    

    S_best.s = (int*)malloc(N * sizeof(int));
    S_best.SizeG = (int*)malloc(M * sizeof(int));
        
    tmpUB = (int*)malloc(M * sizeof(int));
    LBGroup = (int*)malloc(M * sizeof(int));
    UBGroup = (int*)malloc(M * sizeof(int));
    BigThanLB = (int*)malloc(M * sizeof(int));
    vectorElement = (int*)malloc(N * sizeof(int));
    groupElement = (int*)malloc(M * sizeof(int));
    SelectEle = (int*)malloc(N * sizeof(int));
    SelectGroup = (int*)malloc(M * sizeof(int));
    tmpEle = (int*)malloc(N * sizeof(int));
    s1 = (int*)malloc(N * sizeof(int));
    s2 = (int*)malloc(N * sizeof(int));

    p = (int*)malloc(N * sizeof(int));
}

void FIFR_ReleaseMemoryDiversity(void) {
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

    free(p); p = NULL;

}

int FIFR_random_int(int max) {
  GetRNGstate();
  double my_number = unif_rand();
  PutRNGstate();
  return (int) floor(my_number * max);
}
