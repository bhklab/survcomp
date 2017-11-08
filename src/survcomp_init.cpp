#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern "C" {

    void concordanceIndexC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

    /* .Call calls */
    SEXP get_concordanceIndex_onevariable(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

    SEXP mrmr_cIndex(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

    SEXP mrmr_cIndex_ensemble_remove(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


    static const R_CMethodDef CEntries[] = {
        {"concordanceIndexC", (DL_FUNC) &concordanceIndexC, 16},
        {NULL, NULL, 0}
    };


    static const R_CallMethodDef CallEntries[] = {
        {"get_concordanceIndex_onevariable", (DL_FUNC) &get_concordanceIndex_onevariable, 12},
        {"mrmr_cIndex",                      (DL_FUNC) &mrmr_cIndex,                       6},
        {"mrmr_cIndex_ensemble_remove",      (DL_FUNC) &mrmr_cIndex_ensemble_remove,      21},
        {NULL, NULL, 0}
    };

    void R_init_survcomp(DllInfo *dll)
    {
        R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);
        R_forceSymbols(dll, TRUE);
    }

}
