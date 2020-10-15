// Declare Functions
#pragma once

void bicriterion_iterated_local_search_call(
        double *distances, 
        int *N, 
        int *R,
        int *upper_bound, 
        int *WL, 
        double *W, 
        double *Xi, 
        int *partition,
        double *result
);
struct Pareto_element* multistart_bicriterion_pairwise_interchange(
        size_t N, 
        double matrix[N][N], 
        int R, 
        int WL, 
        double weights[WL], 
        int *partition
);

struct Pareto_element* bicriterion_iterated_local_search(struct Pareto_element* head, size_t N, double matrix[N][N], int G, int WL, double weights[WL], double neighbor_percent[2]);
double sample(size_t array_size,double array[array_size]);
void random_partition(size_t N, int G, int partition[N]);
double get_diversity(size_t N, int partition[N], double matrix[N][N]);
double get_dispersion(size_t N, int partition[N], double matrix[N][N]);
void cluster_swap(size_t N, int i, int j, int partition[N]);
void update_pareto(struct Pareto_element** head_ref, size_t N, int partition[N], double diversity, double dispersion);
void compress(size_t N, int partition[N]);
bool paretodominated(struct Pareto_element* head, double diversity, double dispersion);
bool paretoincluded(struct Pareto_element* head, size_t N, int partition[N]);
void push(struct Pareto_element** head_ref, double diversity, double dispersion, size_t N, int partition[N]);
void delete_outdated(struct Pareto_element** head_ref, double diversity, double dispersion);
void linked_list_sample(struct Pareto_element* head, size_t N, int* partition);
int linked_list_length(struct Pareto_element* head);
double random_in_range(double min, double max);
double get_diversity_fast(double diversity, int x,int y, size_t N, int partition[N], double matrix[N][N]);
double get_dispersion_fast(double dispersion, int x,int y, size_t N, int partition[N], double matrix[N][N]);
