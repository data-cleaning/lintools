
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "sparseConstraints.h"
#include "spa.h"
#include "sc_arith.h"
#include "maxdist.h"



static void update_x_k(SparseConstraints *E, double *x, double *w, double *wa, double *alpha, double awa, int k, double *conv){
   
   double *ak = E->A[k];
   int *I = E->index[k];
   int nrag = E->nrag[k];
   
   double ax = 0;
   for ( int j=0; j<nrag; j++ ){
      ax += ak[j] * x[I[j]];
      wa[j] = w[I[j]] * ak[j];
   }

   conv[k] = (ax - E->b[k])/awa;

   double fact = conv[k];
   if ( k >= E->neq ){ // are we at an inequation?
      double alpha_old = alpha[k];
      alpha[k] +=  conv[k];
      if ( alpha[k] < 0 ){
         alpha[k] = 0;
      }
      fact = alpha[k] - alpha_old;
   }
   
   for( int j=0; j < nrag; j++ ){
      x[I[j]] -= wa[j]*fact;
   }
   
}

static void set_zero(double *x, int n){
   for ( int i=0; i<n; ++i ){ 
      x[i] = 0.0;
   }
}

/* Successive projection algorithm, notes.
 *
 * Minimizes x in (x-x0)'W(x-x0) such that Ax <= b holds.
 *
 * exit status: 
 * 0 : ok
 * 1 : not enough memory
 * 2 : divergence
 * 3 : max iterations exceeded
 *
 *  NOTE: C99 std. mentions that calloc's initialization to all bits zero need
 *  not imply a zero double representation on all platforms. We therefore
 *  initialize all (double *)'s ourselves 
 *  
 *  */
int solve_sc_spa(SparseConstraints *E, double *w, double *tol, int *maxiter, double *x  ){
  
   int m = E->nconstraints;
   int n = E->nvar;

   int nrag;
   int niter = 0;
   double *awa    = (double *) malloc(m * sizeof(double));
   double *xw     = (double *) malloc(n * sizeof(double));
   double *alpha  = (double *) malloc(m * sizeof(double));
   double *conv   = (double *) malloc(m * sizeof(double));
   int maxrag = get_max_nrag(E);
   double *wa     = (double *) malloc(maxrag * sizeof(double));

   if ( awa == NULL ||  xw == NULL || alpha == NULL || conv == NULL || wa == NULL ){ 
      // cleanup if one of the objects could nog be allocated
      free(awa); 
      free(xw); 
      free(alpha); 
      free(conv); 
      free(wa);
      return 1;
   } else {
      set_zero(awa,m);
      set_zero(xw,n);
      set_zero(alpha,m);
      set_zero(conv,m);
      set_zero(wa,maxrag);
   }

   int exit_status = 0;
   // we only need w's inverse.
   for ( int k=0; k < n; ++k ){ 
      xw[k] = 1.0/w[k];
   }
   // determine inner products A'W^(-1)A
   for ( int k=0; k < m; k++){
      awa[k] = 0;
      nrag = E->nrag[k];
      for ( int j=0; j<nrag; j++){
         awa[k] += E->A[k][j] * xw[E->index[k][j]] * E->A[k][j];
      }
   }

   // Iterate until convergence, max iterations or divergence detection.
   double diff=DBL_MAX;
   while ( diff > *tol && niter < *maxiter ){

      for ( int k=0; k<m; k++ ) update_x_k(E, x, xw, wa, alpha, awa[k], k, conv);
      ++niter;

      if ( diverged(x,n) || diverged(alpha,m) ){
         exit_status = 2; 
         break;
      }
      // compute convergence criterion
      diff = absmax(conv, awa, E->neq, E->nconstraints); 

   }
   // number of iterations exceeded without convergence?
   if (exit_status != 2 && niter == *maxiter && diff > *tol ) exit_status = 3;

   *tol = sc_diffmax(E,x); // actual difference in current vector
   *maxiter = niter;
   free(wa); 
   free(awa); 
   free(xw); 
   free(alpha); 
   free(conv);
   return exit_status;
}





