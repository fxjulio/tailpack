#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
/*
 * 
 * name: solve
 * @param SEXP A
 * @param double *x
 * 
 * @return x <- A^(-1)*x
 * 
 */

void solve( SEXP A, double *x ){
  int n = Rf_nrows(A);
  int *ipiv = (int *) Calloc(n, int);
  int one = 1;
  char *uplo = "L";
  
  double *avals = (double *) Calloc( n * n, double);
  Memcpy(avals, REAL(A), (size_t)n * n);
  int info, lwork = n+1;
  double *work = (double *) Calloc( lwork, double);
 
  F77_CALL(dsysv)(uplo, &n, &one, avals, &n, ipiv, x, &n, work, &lwork, &info);
  
  Free(ipiv);
  Free(avals);
  Free(work);
  
  if( info < 0 )
    error("argument %d of Lapack routine %s had invalid value",
        -info, "dsysv");
  if (info > 0){
	  Rprintf("A[0,0]:%f A[1,1]:%f\n", REAL(A)[0], REAL(A)[n]);
    error("Lapack routine %s: system is exactly singular: U[%d,%d] = 0",
      "dsysv", info, info);
     
  }
      
}
