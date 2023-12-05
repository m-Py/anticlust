#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* generate this file via
 * tools::package_native_routine_registration_skeleton(".")
 */

/* .C calls */
extern void bicriterion_iterated_local_search_call(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dispersion_anticlustering(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void distance_anticlustering(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kmeans_anticlustering(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fast_kmeans_anticlustering(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"bicriterion_iterated_local_search_call", (DL_FUNC) &bicriterion_iterated_local_search_call,  14},
  {"dispersion_anticlustering",              (DL_FUNC) &dispersion_anticlustering,               9},
  {"distance_anticlustering",                (DL_FUNC) &distance_anticlustering,                 14},
  {"kmeans_anticlustering",                  (DL_FUNC) &kmeans_anticlustering,                  11},
  {"fast_kmeans_anticlustering",                  (DL_FUNC) &fast_kmeans_anticlustering,         8},
  {NULL, NULL, 0}
};

void R_init_anticlust(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
