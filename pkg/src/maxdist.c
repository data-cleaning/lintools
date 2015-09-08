

#include <math.h>


double absmax(double *conv, double *awa, int neq, int nconstraints){

   double d, dmax=0;
   for ( int i=0;  i<neq; ++i ){
      d = fabs(conv[i] * awa[i]);
      if ( d > dmax ) dmax = d;
   }
   for ( int i=neq; i<nconstraints; ++i){
      d = (conv[i] < 0) ? 0 : conv[i]*awa[i];
      if ( d > dmax ) dmax = d;
   }
   return dmax;
}

// d = ||Ax-b||_Inf 
double dc_diffmax(double *A, double *b, double *x, int neq, int m, int n){
   double ax, d, dmax=0;
   for ( int i=0; i < m; ++i ){
      ax = 0;
      for ( int j=0; j < n; ++j ){
         ax += A[i + m*j]*x[j];
      }
      d = ax - b[i];
      if ( i < neq ) d = fabs(d);
      if ( d > dmax ) dmax = d;
   }
   // if all equalities and inequalities are fulfilled
   if ( dmax < 0 ) dmax = 0;
   return dmax;
}


int diverged(double *x, int n){
   for (int i=0; i<n; ++i){
      if (!isfinite(x[i])) return 1;
   }
   return 0;
}




