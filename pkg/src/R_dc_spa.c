
#include <R.h>
#include <Rdefines.h>
#include "dc_spa.h"

SEXP R_dc_solve(SEXP A, SEXP b, SEXP w, SEXP neq, SEXP tol, SEXP maxiter, SEXP x){
   

   SEXP dim;
   double *xx = REAL(x);

   PROTECT(dim = getAttrib(A, R_DimSymbol));

   int m = INTEGER(dim)[0];
   int n = INTEGER(dim)[1];

   if ( n != length(x)) error("%s\n","Number of columns in constraint matrix does not match dimension of x.");
   if ( m != length(b)) error("%s\n","Number of rows in constraint matrix does not mathch dimension of b.");

   // make copies to avoid writing in user space
   SEXP tx;
   PROTECT(tx = allocVector(REALSXP, n));
   for ( int j=0; j<n; j++) REAL(tx)[j] = xx[j];

   double xtol = REAL(tol)[0];
   int xmaxiter = INTEGER(maxiter)[0];
   

   int s = dc_solve(
      REAL(A), 
      REAL(b), 
      REAL(w), 
      m, 
      n, 
      INTEGER(neq)[0],
      &xtol,
      &xmaxiter,
      REAL(tx)
   );
   
   SEXP status, niter, eps;
   PROTECT(status= allocVector(INTSXP,1));
   PROTECT(niter = allocVector(INTSXP,1));
   PROTECT(eps   = allocVector(REALSXP,1));

   INTEGER(status)[0] = s;
   INTEGER(niter)[0]  = xmaxiter;
   REAL(eps)[0] = xtol;
   

   setAttrib(tx, install("status"),status);
   setAttrib(tx, install("niter"), niter);
   setAttrib(tx, install("tol"), eps);

   UNPROTECT(5);
   return tx;
}



