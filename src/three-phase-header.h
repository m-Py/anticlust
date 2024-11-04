#ifndef THREE_PHASE_HEADER_H 
#define THREE_PHASE_HEADER_H

typedef struct {
    int *s;        // cluster membership of each element
    int *SizeG;    // size of each cluster
    double objective;   // global cost function value of solution
} Solution;

int random_int(int max);
double uniform_rnd_number(void);
void swap_elements(int* a, int* b);
void fisher_yates_shuffle(int arr[], int n);
void RandomInitialSol(int s[], int SizeG[]);
void UndirectedPerturbation(int theta, int partition[], int SizeGroup[]);
double LocalSearchCriterionCalculation(Solution* sol1, Solution* sol2);
int CompareSolution(const void *first, const void *second);

#endif // THREE_PHASE_HEADER_H