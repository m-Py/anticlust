#ifndef FIFR_HEADER_H
#define FIFR_HEADER_H


typedef struct Solution {
    int* s;
    int* SizeG;
    double objective;
} Solution;

int  FIFR_random_int(int max);
void FIFR_swap_elements(int *a, int *b);
void FIFR_fisher_yates_shuffle(int arr[], int n);

int  FIFR_CompareSolution(const void *first, const void *second);

double FIFR_LocalSearchCriterionCalculation(Solution* sol1, Solution* sol2);

void FIFR_ClearDeltaMatrix(void);
void FIFR_BuildDeltaMatrix(int partition[]);
void FIFR_OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup);

void FIFR_BuildGroupDiversityForCrossover(int partition[], double groupDiversity[]);
void FIFR_process_partition(
    double* groupDiversity,
    int* partition,
    int* tmpUB,
    int* childSolution,
    int* vectorElement,
    int M,
    int N,
    int* elementCount,
    int* targetGroup
);
void FIFR_CrossoverDiversity(int partition1[], int partition2[], int childSolution[], int scSizeGroup[]);

void FIFR_SearchAlgorithmDiversity(void);

double uniform_rnd_number(void);
void InitialSolution(int s[], int SizeG[]);
void Perturbation(int L, int s[], int SizeGroup[]);

#endif //FIFR_HEADER_H