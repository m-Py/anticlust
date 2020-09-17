#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void distance_anticlustering(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kmeans_anticlustering(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"distance_anticlustering", (DL_FUNC) &distance_anticlustering,  9},
  {"kmeans_anticlustering",   (DL_FUNC) &kmeans_anticlustering,   11},
  {NULL, NULL, 0}
};

void R_init_anticlust(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
