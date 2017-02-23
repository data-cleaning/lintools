#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP all_finite_double(SEXP);
extern SEXP R_dc_solve(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_get_nconstraints(SEXP);
extern SEXP R_get_nvar(SEXP);
extern SEXP R_print_sc(SEXP, SEXP, SEXP);
extern SEXP R_sc_diffmax(SEXP, SEXP);
extern SEXP R_sc_diffsum(SEXP, SEXP);
extern SEXP R_sc_diffvec(SEXP, SEXP);
extern SEXP R_sc_from_sparse_matrix(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_sc_multvec(SEXP, SEXP);
extern SEXP R_solve_sc_spa(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"all_finite_double",       (DL_FUNC) &all_finite_double,       1},
    {"R_dc_solve",              (DL_FUNC) &R_dc_solve,              7},
    {"R_get_nconstraints",      (DL_FUNC) &R_get_nconstraints,      1},
    {"R_get_nvar",              (DL_FUNC) &R_get_nvar,              1},
    {"R_print_sc",              (DL_FUNC) &R_print_sc,              3},
    {"R_sc_diffmax",            (DL_FUNC) &R_sc_diffmax,            2},
    {"R_sc_diffsum",            (DL_FUNC) &R_sc_diffsum,            2},
    {"R_sc_diffvec",            (DL_FUNC) &R_sc_diffvec,            2},
    {"R_sc_from_sparse_matrix", (DL_FUNC) &R_sc_from_sparse_matrix, 5},
    {"R_sc_multvec",            (DL_FUNC) &R_sc_multvec,            2},
    {"R_solve_sc_spa",          (DL_FUNC) &R_solve_sc_spa,          5},
    {NULL, NULL, 0}
};

void R_init_lintools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
