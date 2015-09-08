
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "sparseConstraints.h"

static void set_null_dbl(double **x, int n){
   for ( int i=0; i<n; ++i ){
      x[i] = NULL;
   }
}

static void set_null_int(int **x, int n){
   for ( int i=0; i<n; ++i ){
      x[i] = NULL;
   }
}

SparseConstraints * sc_new( int m ){
    
   SparseConstraints *E;
   E  = (SparseConstraints *) calloc(1, sizeof(SparseConstraints));
   if ( E == NULL ){ 
      return NULL;
   }
   E->nconstraints = m;
   E->A     = (double **) calloc(E->nconstraints, sizeof(double *));
   E->index = (int **) calloc(E->nconstraints, sizeof(int *));
   E->nrag  = (int *) calloc(E->nconstraints, sizeof(int));
   E->b     = (double *) calloc(E->nconstraints, sizeof(double));

   if ( E->A == NULL || E->index == NULL || E->nrag == NULL || E->b == NULL ){
      free(E->A);
      free(E->index);
      free(E->nrag);
      free(E->b);
      return NULL;
   } else {
      set_null_dbl(E->A, E->nconstraints);
      set_null_int(E->index, E->nconstraints);
      return E;
   }
}

void sc_del(SparseConstraints *E){

   if ( E == NULL ) return;
   for ( int i = 0; i < E->nconstraints; i++){
      free(E->A[i]);
      free(E->index[i]);
   }
   free(E->b);
   free(E->nrag);
   free(E->A);
   free(E->index);
   free(E);
}


static int get_row_end(int *rows, int nrows, int row_start){
    int row_nr = rows[row_start];
    int row_end = row_start + 1;
    while (row_end < (nrows) && rows[row_end] == row_nr) row_end++;
    return row_end;
}



/* Generates a sparse representation of Ax <op> b, where the first neq <op> are
 * '=' and the rest are interpreted as '<='. 
 *
 * - A is specified in row, column, coefficient triples.
 * - no dimension checking on A or b is performed.
 * - It is assumed that  int *rows is sorted in increasing order!
 */
SparseConstraints * sc_from_sparse_matrix(int *rows, int *cols, double *coef, int ncoef, double *b, int m, int neq ){

   int maxcol=0, k;
   int row_start = 0, row_end; 

   SparseConstraints *E = sc_new(m);

   if ( E == NULL ) return NULL;

   

   for ( int irow=0; irow < m; irow++){
      E->b[irow] = b[irow];
      row_end = get_row_end(rows, ncoef, row_start);
     
      E->nrag[irow]    = row_end - row_start;
      E->index[irow]   = (int *) calloc( E->nrag[irow], sizeof(int));
      E->A[irow]       = (double *) calloc( E->nrag[irow], sizeof(double));
      if ( E->A[irow] == NULL || E->index[irow] == NULL ){ // no memory = no joy 
         sc_del(E);
         return NULL;
      }            

      k = 0;
      for ( int j=row_start; j < row_end; j++){
         E->A[irow][k] = coef[j];
         E->index[irow][k] = cols[j];
         ++k;
         if (cols[j] > maxcol) maxcol = cols[j];
      }
      row_start = row_end;
   } 
   
   E->neq = neq;
   E->nvar = maxcol+1;

   return E;

}

int get_max_nrag(SparseConstraints *E){
   int nmax = INT_MIN;
   for ( int i=0; i < E->nconstraints; ++i ){
      if ( nmax < E->nrag[i] ) nmax = E->nrag[i];
   }
  return nmax;
}




