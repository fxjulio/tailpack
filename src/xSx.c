#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

//S simetrica [nxn]
double xSx( SEXP S, double *x  ){
  int n = Rf_nrows(S);
  int i, j;
  double ans = 0.0, *ps = REAL(S);

  for( i=0; i<n; i++ ){
    ans += ps[ (n+1)*i ]*x[i]*x[i];
    for( j=0; j<i; j++ ){
      ans += 2.0*ps[ n*j + i ]*x[i]*x[j];
    }
    
  }

  return ans;
}
