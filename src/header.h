// Declare Functions
#pragma once

void bicriterion_iterated_local_search_call(
        double *distances, 
        double *disp_distances,
        int *N, 
        int *R,
        int *upper_bound, 
        int *WL, 
        double *W, 
        double *Xi, 
        int *partition,
        int *frequencies, // frequency of each partition
        int *use_init_partitions,
        int *init_partitions,
        int *result,
        int *mem_error
);
struct Pareto_element* multistart_bicriterion_pairwise_interchange(
        size_t N, 
        double matrix[N][N], 
        double matrix2[N][N], 
        size_t R, 
        size_t WL, 
        double weights[WL], 
        int *partition, int *frequencies, int *use_init_partitions,
        int *init_partitions
);

void shuffle_permutation(int N, int *permutation);

struct Pareto_element* bicriterion_iterated_local_search(
        struct Pareto_element* head, 
        size_t N, 
        double matrix[N][N], 
        double matrix2[N][N], 
        size_t G, 
        size_t WL, 
        double weights[WL], 
        double neighbor_percent[2],
        int *frequencies
);
double sample(size_t array_size,double array[array_size]);
double get_diversity(size_t N, int* partition, double matrix[N][N], int *frequencies);
double get_dispersion(size_t N, int* partition, double matrix[N][N]);
void cluster_swap(size_t i, size_t j, int* partition);
int update_pareto(struct Pareto_element** head_ref, size_t N, int* partition, double diversity, double dispersion);
bool paretodominated(struct Pareto_element* head, double diversity, double dispersion);
int push(struct Pareto_element** head_ref, double diversity, double dispersion, size_t N, int* partition);
void delete_outdated(struct Pareto_element** head_ref, double diversity, double dispersion);
void linked_list_sample(struct Pareto_element* head, size_t N, int* partition);
int linked_list_length(struct Pareto_element* head);
double random_in_range(double min, double max);
double get_diversity_fast(double diversity, int x,int y, size_t N, int* partition, double matrix[N][N], int *frequencies);
double get_dispersion_fast(double dispersion, int x,int y, size_t N, int* partition, double matrix[N][N]);
void free_pareto_set(struct Pareto_element* head);
double uniform_rnd_number(void);
double uni_rnd_number_range(double min, double max);
int random_integer(int min, int max);
