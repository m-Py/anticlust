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
struct node* fill_list(
        struct node *HEAD, 
        struct element *data,
        size_t i
);
void print_cluster(
        struct node *cluster, 
        size_t M
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
double update_objective(
        size_t k, 
        size_t m, 
        double centers[k][m], 
        size_t cl1, 
        size_t cl2, 
        struct node *head_array[k], 
        double VAR_OBJECTIVE[k]
);
void print_elements(
        size_t n, 
        size_t m, 
        struct node **ptr_to_nodes
);
void cp_array(
        size_t n, 
        double origin[n], 
        double target[n]
);
void cp_matrix(
        size_t n, 
        size_t m, 
        double origin[n][m],
        double target[n][m]
);
