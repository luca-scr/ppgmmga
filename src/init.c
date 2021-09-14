#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
Routines registration obtained with 

tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
 
FIXME: Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ppgmmga_encodebasis(SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_EntropyGauss(SEXP, SEXP);
extern SEXP _ppgmmga_EntropyGMM(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_EntropyMCapprox(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_EntropySOTE(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_EntropyUT(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_EntropyVAR(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_LinTransf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ppgmmga_logsumexp(SEXP);
extern SEXP _ppgmmga_orth(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ppgmmga_encodebasis",     (DL_FUNC) &_ppgmmga_encodebasis,     3},
    {"_ppgmmga_EntropyGauss",    (DL_FUNC) &_ppgmmga_EntropyGauss,    2},
    {"_ppgmmga_EntropyGMM",      (DL_FUNC) &_ppgmmga_EntropyGMM,      4},
    {"_ppgmmga_EntropyMCapprox", (DL_FUNC) &_ppgmmga_EntropyMCapprox, 5},
    {"_ppgmmga_EntropySOTE",     (DL_FUNC) &_ppgmmga_EntropySOTE,     5},
    {"_ppgmmga_EntropyUT",       (DL_FUNC) &_ppgmmga_EntropyUT,       5},
    {"_ppgmmga_EntropyVAR",      (DL_FUNC) &_ppgmmga_EntropyVAR,      5},
    {"_ppgmmga_LinTransf",       (DL_FUNC) &_ppgmmga_LinTransf,       6},
    {"_ppgmmga_logsumexp",       (DL_FUNC) &_ppgmmga_logsumexp,       1},
    {"_ppgmmga_orth",            (DL_FUNC) &_ppgmmga_orth,            2},
    {NULL, NULL, 0}
};

void R_init_ppgmmga(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
