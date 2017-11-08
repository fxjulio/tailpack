#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
/*
 * 
 * name: rmvnorm
 * @param mu
 * @param P
 * @param x
 * @param jobz
 * 
 * @return 
 * 
 */

double rmvnorm( double *mu, SEXP P, double *x, char *jobz ){
  char *uplo = "L", *trans = "N";
  int inc = 1;
  int p = Rf_nrows(P);
  double *eigvec = (double *) Calloc( p * p, double);
  Memcpy(eigvec, REAL(P), (size_t)p * p);
  double *eigval = (double *) Calloc( p, double);
  int lwork = 3*p-1;
  double *work = (double *) Calloc( lwork, double);
  int info;
  F77_CALL(dsyev)( jobz, uplo, &p, eigvec, &p, eigval, work, &lwork, &info);
  
  if( info != 0 ){
    Free(eigvec);
    Free(eigval);
    Free(work);
  }
  
  if( info < 0 )
    error("argument %d of Lapack routine %s had invalid value",
        -info, "dsyev");
  if (info > 0)
    error("Lapack routine %s: the algorithm failed to converge; %d off-diagonal elements of an intermediate tridiagonal form did not converge to zero.",
      "dsyev", info);
      
  if( eigval[0] < 0.0 ){
    Rprintf("Matrix is not positive-semidefinite eigval[0]: %f", eigval[0]);
    return eigval[0];
  }
  
  if( strcmp(jobz, "N") == 0){
    double ldetP = 0.0;
    for( int i=0; i<p; i++ ){
      ldetP += log(eigval[i]);
    }
    Free(eigvec);
    Free(eigval);
    Free(work);
    return ldetP;
  }  
  double sum;
  double *R = (double *) Calloc( p * p, double); // R=P^(-1/2)P^(-1/2)
  for( int i=0; i < p ; i++ ){
    sum = 0.0;
    for( int k=0; k < p; k++)
      sum += 1.0/sqrt(eigval[k])*eigvec[ p*k + i]*eigvec[ p*k + i] ;
    R[ (p+1)*i ] = sum;
    for( int j=0; j < i ; j++){
      sum = 0.0;
      for( int k=0; k < p; k++)
        sum += 1.0/sqrt(eigval[k])*eigvec[ p*k +i]*eigvec[ p*k +j] ;
      R[ p*i + j ] = sum;
      R[ p*j + i ] = sum;
    }
    
    
  }

  double *Z = (double *) Calloc( p, double);
  
  double ldetP = 0.0;
  for( int i=0; i<p; i++ ){
    Z[i] = rnorm(0,1);
    x[i] = mu[i];
    ldetP += log(eigval[i]);
  }
  
  double blas_alpha = 1.0;
  double blas_beta = 1.0;
  F77_NAME(dgemv)(trans, &p, &p, &blas_alpha, R, &p, Z, &inc,
    &blas_beta, x, &inc);
    
  Free(eigvec);
  Free(eigval);
  Free(work);
  Free(R);
  Free(Z);
  return ldetP;
}
