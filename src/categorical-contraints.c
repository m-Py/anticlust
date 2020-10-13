 
#include "declarations.h"
#include <stdlib.h> 
 
/* Determine number of categories 
 * Is either 1 - if argument *USE_CATS is FALSE - or *C
 */
size_t number_of_categories(int *USE_CATS, int *C) {
        if (*USE_CATS) {
                return (size_t) *C;
        }
        return 1;
}

/* For each category, get the number of members in it
 * If no categories are used (i.e., *USE_CATS = FALSE), the program
 * proceed as if there is one big category that contains every 
 * element, i.e., n elements.
 */
int get_cat_frequencies(int *USE_CATS, int *CAT_frequencies, size_t n) {
        if (*USE_CATS) {
                return *CAT_frequencies;
        }
        return n;
}

/* Function that writes to a "global" variable `CATEGORY_HEADS`
 *   
 *   `CATEGORY_HEADS` is an array of pointers to arrays (each entry points 
 *   to an array). The i'th element of `C_HEADS` contains an array of all 
 *   indices of the elements that have the i'th category.
 *
 * If no categorical constraints are induced (that is, *USE_CATS = FALSE), 
 * `CATEGORY_HEADS` points to a single array that contains all indices 1, ..., n.
 * This indicates that all elements serve as exchange partners in the exchange
 * algorithm (in general, only elements within the same array of `CATEGORY_HEADS`
 * serve as exchange partners).
 * 
 */
int get_indices_by_category(size_t n, size_t c, size_t *CATEGORY_HEADS[c], 
                            int *USE_CATS, 
                            int *categories, int *CAT_frequencies, 
                            struct element POINTS[n]) {
        // Set up array of exchange partners
        if (*USE_CATS) {
                if (set_up_categories_list(n, c, POINTS, CATEGORY_HEADS, 
                                     categories, CAT_frequencies) == 1) {
                        return 1;
                }
                return 0;
        }
        CATEGORY_HEADS[0] = (size_t*) malloc(n * sizeof(size_t));
        if (CATEGORY_HEADS[0] == NULL) {
                return 1;
        }
        for (size_t i = 0; i < n; i++) {
                CATEGORY_HEADS[0][i] = i;
        }
        return 0;
}

/* This function actually does the work of creating the index arrays pointed
 * to by `C_HEADS`. Sets up a cluster structure where each category corresponds
 * to a cluster. Then, proceeds to read out the indices of all elements by category. 
 */
int set_up_categories_list(size_t n, size_t c, struct element POINTS[n], 
                     size_t *CATEGORY_HEADS[c], int *categories, 
                     int *CAT_frequencies) {
        
        struct node *HEADS[c]; // used for filling `CATEGORY_HEADS` 
        if (initialize_cluster_heads(c, HEADS) == 1) {
                return 1; 
        }
        
        // Set up array of pointers-to-nodes, return if memory runs out
        struct node *PTR_NODES[n];
        if (fill_cluster_lists(n, c, categories, POINTS, PTR_NODES, HEADS) == 1) {
                free_cluster_list(c, HEADS, c);
                return 1;
        }
        
        // Initialize the index arrays
        struct node *tmp; 
        for (size_t i = 0; i < c; i++) {
                size_t n_cats = (size_t) CAT_frequencies[i];
                CATEGORY_HEADS[i] = (size_t*) malloc(n_cats * sizeof(size_t));
                if (CATEGORY_HEADS[i] == NULL) {
                        free_category_indices(c, CATEGORY_HEADS, i);
                        free_cluster_list(c, HEADS, c);
                        return 1;
                }
                // Now write `CATEGORY_HEADS`! Fills all `c` arrays with indices, 
                // based on the category lists `HEAD[0], HEAD[1], ..., HEAD[c]`
                tmp = HEADS[i]->next;
                size_t j = 0;
                while (tmp != NULL) {
                        CATEGORY_HEADS[i][j] = tmp->data->ID;
                        j++;
                        tmp = tmp->next;
                } 
        }

        // free temporary category lists
        free_cluster_list(c, HEADS, c);
        return 0;
}
