
#include <R.h>
#include <Rdefines.h>
#include "sparseConstraints.h"
#include "sc_arith.h"


SEXP R_sc_multvec(SEXP p, SEXP x){
   PROTECT(p);
   PROTECT(x);

   SparseConstraints *xp = R_ExternalPtrAddr(p);

   SEXP Ax; 
   PROTECT(Ax = allocVector(REALSXP, xp->nconstraints));
   
   sc_multvec(xp, REAL(x), REAL(Ax));

   UNPROTECT(3);
   return Ax;

}


SEXP R_sc_diffvec(SEXP p, SEXP x){
   PROTECT(p);
   PROTECT(x);
   SparseConstraints *xp = R_ExternalPtrAddr(p);

   SEXP dv;
   PROTECT(dv = allocVector(REALSXP, xp->nconstraints));

   sc_diffvec(xp, REAL(x), REAL(dv));

   UNPROTECT(3);
   return dv;

}

SEXP R_sc_diffmax(SEXP p, SEXP x){

   PROTECT(p);
   PROTECT(x);
   SparseConstraints *xp = R_ExternalPtrAddr(p);

   SEXP d;
   PROTECT(d = allocVector(REALSXP,1));

   REAL(d)[0] = sc_diffmax( xp, REAL(x) );
   UNPROTECT(3);

   return d;

}

SEXP R_sc_diffsum(SEXP p, SEXP x){

   PROTECT(p);
   PROTECT(x);
   SparseConstraints *xp = R_ExternalPtrAddr(p);

   SEXP d;
   PROTECT(d = allocVector(REALSXP,1));

   REAL(d)[0] = sc_diffsum( xp, REAL(x) );
   UNPROTECT(3);

   return d;

}

