
#include <R.h>
#include <Rdefines.h>
#include "sparseConstraints.h"
#include "spa.h"


SEXP R_solve_sc_spa(SEXP p, SEXP x, SEXP w, SEXP tol, SEXP maxiter){

   PROTECT(p);
   PROTECT(x);
   PROTECT(w);
   PROTECT(tol);
   PROTECT(maxiter);
   SEXP niter, eps, status;
   SparseConstraints *xp = R_ExternalPtrAddr(p);
    
   // make copies outside R to prevent writing in userspace.
   double xtol = REAL(tol)[0];
   int xmaxiter = INTEGER(maxiter)[0];
   int s;
   double *xx = REAL(x);
   SEXP tx;

   PROTECT(tx = allocVector(REALSXP, length(x)));
   for ( int i=0; i<length(x); i++) REAL(tx)[i] = xx[i];

   // solve
   s = solve_sc_spa(xp, REAL(w) , &xtol, &xmaxiter, REAL(tx)); 

   // return answer to R
   PROTECT(status = allocVector(INTSXP,1));
   PROTECT(niter = allocVector(INTSXP,1));
   PROTECT(eps = allocVector(REALSXP,1));

   INTEGER(status)[0] = s;
   INTEGER(niter)[0] = xmaxiter;
   REAL(eps)[0] = xtol;

   setAttrib(tx,install("niter"), niter);
   setAttrib(tx,install("tol"), eps);
   setAttrib(tx,install("status"), status);

   UNPROTECT(9);
   return tx;
}



