#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double logistic( double x){
  double ans;
  if( x > 16.0 ) x = 16.0;
  ans = 1.0 / ( 1.0 + exp(-x));
  if( ans > 0.999 ) ans = 0.999;
  if( ans < 0.001 ) ans = 0.001;
  
  return ans;
  
}

//Devuelve una matriz simetrica t(X)WX + beta Y  donde W es diagonal
void tXWX( SEXP X, double *w, double beta, SEXP Y, SEXP ans ){
  int n = Rf_nrows(X);
  int p = Rf_ncols(X);
  int i, j, k;
  double sum, *pans = REAL(ans), *px = REAL(X), *py = REAL(Y);
  
  for( i=0; i<p; i++ ){
    for( j=0; j<p; j++ ){
      sum = beta * py[ p*j + i ];
      for( k=0; k < n ; k++){
        sum += w[k] * px[ n*i + k ] * px[ n*j + k];
      }
      pans[ p*j + i ] = sum;
    }
  }
  
}

/*  */

void tXWytilde( SEXP X, double *w, double *ytilde, double *ans ){
  int n = Rf_nrows(X);
  int p = Rf_ncols(X);
  double sum, *px = REAL(X);
  int i, j;
  
  for( i=0; i < p; i++){
    sum = 0.0;
    for( j=0; j < n; j++){
      sum += px[ n*i + j] * w[j] * ytilde[j];
    }
    ans[i] = sum;
  }
  
}
