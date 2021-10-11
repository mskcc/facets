#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void F77_SUB(scansnp)(int *n, double *maploc, double *het, double *keep, double *nbhd);
void F77_SUB(t2maxo)(int *n, double *sx, int *ihet, int *iseg, double *ostat, int *nhet, double *rnij, double *rhij, double *delij);

static const R_FortranMethodDef FortranEntries[] = {
    {"scansnp", (DL_FUNC) &F77_SUB(scansnp), 5},
    {"t2maxo",  (DL_FUNC) &F77_SUB(t2maxo),  9},
    {NULL, NULL, 0}
};

void R_init_facets(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
