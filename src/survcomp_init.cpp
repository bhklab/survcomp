#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void concordanceIndexC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern "C" SEXP get_concordanceIndex_onevariable(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP mrmr_cIndex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP mrmr_cIndex_ensemble_remove(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"get_concordanceIndex_onevariable", (DL_FUNC) &get_concordanceIndex_onevariable, 12},
    {"mrmr_cIndex",                      (DL_FUNC) &mrmr_cIndex,                       6},
    {"mrmr_cIndex_ensemble_remove",      (DL_FUNC) &mrmr_cIndex_ensemble_remove,      21},
    {NULL, NULL, 0}
};

void R_init_survcomp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
