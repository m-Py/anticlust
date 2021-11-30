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
        int *result,
        int *mem_error
);
struct Pareto_element* multistart_bicriterion_pairwise_interchange(
        size_t N, 
        double matrix[N][N], 
        size_t R, 
        size_t WL, 
        double weights[WL], 
        int *partition
);

void shuffle_permutation(int N, int *permutation);

struct Pareto_element* bicriterion_iterated_local_search(
        struct Pareto_element* head, 
        size_t N, 
        double matrix[N][N], 
        size_t G, 
        size_t WL, 
        double weights[WL], 
        double neighbor_percent[2]
);
double sample(size_t array_size,double array[array_size]);
double get_diversity(size_t N, int* partition, double matrix[N][N]);
double get_dispersion(size_t N, int* partition, double matrix[N][N]);
void cluster_swap(size_t i, size_t j, int* partition);
int update_pareto(struct Pareto_element** head_ref, size_t N, int* partition, double diversity, double dispersion);
bool paretodominated(struct Pareto_element* head, double diversity, double dispersion);
int push(struct Pareto_element** head_ref, double diversity, double dispersion, size_t N, int* partition);
void delete_outdated(struct Pareto_element** head_ref, double diversity, double dispersion);
void linked_list_sample(struct Pareto_element* head, size_t N, int* partition);
int linked_list_length(struct Pareto_element* head);
double random_in_range(double min, double max);
double get_diversity_fast(double diversity, int x,int y, size_t N, int* partition, double matrix[N][N]);
double get_dispersion_fast(double dispersion, int x,int y, size_t N, int* partition, double matrix[N][N]);
void free_pareto_set(struct Pareto_element* head);
double uniform_rnd_number();
double uni_rnd_number_range(double min, double max);
int random_integer(int min, int max);
