

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "maxdist.h"


static void update_x_k(double *A, double *b, double *x, int neq, int m, int n, double *w, double *wa, double *alpha, double awa, int k, double *conv){

   double alpha_old = alpha[k];
   double ax = 0, fact;

   for ( int j=0; j<n; j++){
      ax += A[k + j*m] * x[j];
      wa[j] = w[j] * A[k + j*m];
   }

   conv[k] = (ax - b[k])/awa;
   
   fact = conv[k];
   if ( k >= neq ){
      alpha[k] += conv[k];
      if ( alpha[k] < 0 ) alpha[k] = 0;
      fact = alpha[k] - alpha_old;
   }

   for ( int j=0; j<n; j++ ){
      x[j] -=   wa[j]*fact;
   }

}


// optimal adjustments with dense constraints.
int dc_solve(double *A, double *b, double *w, int m, int n, int neq, double *tol, int *maxiter, double *x){
   
   int niter = 0;

   double *awa = (double *) calloc(m, sizeof(double)); 
   double *xw = (double *) calloc(n, sizeof(double));
   double *alpha = (double *) calloc(m, sizeof(double));
   double *wa = (double *) calloc(n, sizeof(double));
   double *conv = (double *) calloc(m, sizeof(double));

   // in case of emergency: wee haave too get autofhea (Schwartzenegger style)
   if ( awa == NULL || xw == NULL|| alpha == NULL|| conv == NULL || wa == NULL ){ 
      free(awa); 
      free(xw); 
      free(alpha); 
      free(conv); 
      free(wa);
      return 1;
   }
   double diff = DBL_MAX; 
   int exit_status = 0;
   
   // we only need w's inverse.
   for ( int k=0; k < n; ++k ){
      xw[k] = 1/w[k];
   }
   
   // determine diag(A'W^(-1)A)
   for ( int k=0; k < m; k++){
      awa[k] = 0;
      for ( int j=0; j<n; j++ ){
         awa[k] += A[k+j*m]*xw[j]*A[k+j*m];
      }
   }


   while ( diff > *tol && niter < *maxiter ){

      for (int k=0; k<m; k++) update_x_k(A, b, x, neq, m, n, xw, wa, alpha, awa[k], k, conv);
      ++niter;

      if ( diverged(x,n) || diverged(alpha,m) ){
         exit_status = 2; 
         break;
      }
      diff = absmax(conv, awa, neq, m);
   }
   // number of iterations exceeded without convergence?
   if (exit_status != 2 && niter == *maxiter && diff > *tol ){ 
      exit_status = 3;
   }
   *tol = dc_diffmax(A, b, x, neq, m, n); // actual max abs diff.
   *maxiter = niter;

   free(wa); 
   free(awa); 
   free(xw); 
   free(alpha); 
   free(conv);
   return exit_status;
}




