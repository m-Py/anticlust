#pragma once
#include <stdlib.h> 

/* Define struct containing data points */
struct element
{
        size_t ID; // index of the node in original order of the data
        size_t cluster; // index of element's cluster
        double *values; // array of length M, element's data values
        size_t category; // index of element's category
};

/* Define struct for nodes in linked list (representing a cluster) */
struct node
{
        struct element *data; // Pointer to data point
        struct node *next; // pointer to next node
};

// Declare Functions
void kmeans_anticlustering(
        double *data, 
        int *N, 
        int *M, 
        int *K, 
        int *frequencies,
        int *clusters, 
        int *USE_CATS, 
        int *C, 
        int *CAT_frequencies,
        int *categories,
        int *mem_error
);

size_t number_of_categories(int *USE_CATS, int *C);
int get_cat_frequencies(int *USE_CATS, int *CAT_frequencies, size_t n);

double euclidean_squared(
        double *x, 
        double *y, 
        size_t m
);
struct node* append_to_cluster(
        struct node *HEAD, 
        struct element *data
);
void compute_center(
        size_t m, 
        double center[m], 
        struct node *HEAD, 
        int frequency
);
int fill_data_points(
        double *data, 
        size_t n, 
        size_t m, 
        struct element POINTS[n], 
        int *clusters,
        int *USE_CATS,
        int *categories
);
double cluster_var(
        size_t m, 
        struct node *HEAD, 
        double center[m]
);
void swap(
        size_t n, 
        size_t i, 
        size_t j, 
        struct node *PTR_NODES[n]
);
void update_centers(
        size_t k, 
        size_t m, 
        double CENTERS[k][m],
        struct node *one, 
        struct node *two,
        int *frequencies
);
void update_objective_by_cluster(
        size_t k, 
        size_t m, 
        double centers[k][m], 
        size_t cl1, 
        size_t cl2, 
        struct node *HEADS[k], 
        double OBJ_BY_CLUSTER[k]
);
void copy_array(
        size_t n, 
        double origin[n], 
        double target[n]
);
void copy_matrix(
        size_t n, 
        size_t m, 
        double origin[n][m],
        double target[n][m]
);
int initialize_cluster_heads(
        size_t k, 
        struct node *HEADS[k]
);

int fill_cluster_lists(
        size_t n, 
        size_t k,
        int *clusters,
        struct element POINTS[n],
        struct node *PTR_NODES[n], 
        struct node *PTR_CLUSTER_HEADS[k]
);
void objective_by_cluster(
        size_t m, 
        size_t k, 
        double OBJ_BY_CLUSTER[k], 
        double CENTERS[k][m], 
        struct node *HEADS[k]
);
double array_sum(
        size_t k, 
        double ARRAY[k]
);

int get_indices_by_category(
        size_t n, 
        size_t c, 
        size_t *CATEGORY_HEADS[c], 
        int *USE_CATS, 
        int *categories, 
        int *CAT_frequencies, 
        struct element POINTS[n]
);
int set_up_categories_list(
        size_t n, 
        size_t c, 
        struct element POINTS[n], 
        size_t *CATEGORY_HEADS[c], 
        int *categories, 
        int *CAT_frequencies
);

/* Free functions */
void free_points(
        size_t n, 
        struct element POINTS[n],
        size_t i
);

void free_cluster_list(
        size_t k, 
        struct node *PTR_CLUSTER_HEADS[k],
        size_t i
);

void free_category_indices(
        size_t c, 
        size_t *CATEGORY_HEADS[c], 
        size_t i
);

void free_distances(
        size_t n, 
        double *DISTANCES[n], 
        size_t i
);

// For distance anticlustering, objective functions
void distance_objective(
        size_t n, 
        size_t k, 
        double *distances[n], 
        double OBJ_BY_CLUSTER[k], 
        struct node *HEADS[k]
);

double distances_within(
        size_t n, 
        double *distances[n], 
        struct node *HEAD
);

double distances_one_element(
        size_t n, 
        double *distances[n], 
        struct node *start_node, 
        size_t ID
);

// For distance anticlustering, objective functions
double dispersion_objective(
        size_t n, 
        size_t k, 
        double *distances[n], 
        struct node *HEADS[k]
);

double minimun_distance_cluster(
        size_t n, 
        double *distances[n], 
        struct node *HEAD
);

double minimin_distance_element(
        size_t n, 
        double *distances[n], 
        struct node *start_node, 
        size_t ID
);

int has_node_dispersion(
        size_t n, 
        double *distances[n], 
        double dispersion, 
        struct node *HEAD, 
        struct node *node,
        int after
);

// Declare Functions
void fast_kmeans_anticlustering(
        double *data,
        int *N,
        int *M,
        int *K,
        int *frequencies,
        int *clusters,
        int *partners,
        int *k_neighbours
);

void fast_update_one_center(size_t index_removed_from_cluster, 
                            size_t index_added_to_cluster, 
                            size_t n, size_t m, double *data, 
                            double CENTER[m], int frequency);
// for fast new k-means implementation
void init_centers(size_t k, size_t m, size_t n, double CENTERS[k][m], int* clusters, int* frequencies, double* data);
void init_overall_centroid(size_t m, size_t n, double OVERALL_CENTROID[m], double* data);
void fast_swap(int *clusters, size_t i, size_t j);
size_t one_dim_index(size_t i, size_t j, size_t n);
double weighted_array_sum(size_t k, int* frequencies, double ARRAY[k]);
