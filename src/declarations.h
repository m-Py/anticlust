/* Define struct containing data points */
struct element
{
        size_t ID; // index of the node in original order of the data
        size_t cluster; // index of element's cluster
        double *values; // array of length M, element's data values
};

/* Define struct for nodes in linked circular list (representing a cluster) */
struct node
{
        struct element *data; // Pointer to data point
        struct node *next; // pointer to next node
};

// Declare Functions
void anticlustering(
        double *data, 
        int *N, 
        int *M, 
        int *K, 
        int *frequencies, 
        int *clusters
);
double euclidean_squared(
        double *x, 
        double *y, 
        size_t m
);
struct node* append_element_to_cluster(
        struct node *HEAD, 
        struct element *data,
        size_t i
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
        struct element points[n], 
        int *clusters
);
double cluster_variance(
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
        struct node *PTR_CLUSTER_HEADS[k], 
        double VAR_OBJECTIVE[k]
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
        size_t n, 
        size_t k, 
        struct node *PTR_CLUSTER_HEADS[k],
        struct element POINTS[n]
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
        double VAR_OBJECTIVE[k], 
        double CENTERS[k][m], 
        struct node *PTR_CLUSTER_HEADS[k]
);
double array_sum(
        size_t k, 
        double ARRAY[k]
);


/* Free functions */
void free_points(
        size_t n, 
        struct element POINTS[n]
);
void free_nodes(
        size_t k, 
        struct node *PTR_CLUSTER_HEADS[k]
);

void print_memory_error();
