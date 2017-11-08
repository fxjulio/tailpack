#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

double xSx( SEXP , double *  );
void tXWX( SEXP , double *, double, SEXP, SEXP  );
void tXWytilde( SEXP, double *, double *, double * );
double logistic( double );
void solve( SEXP A, double * );
double rmvnorm( double *, SEXP, double *, char *);

SEXP bsmoothExp( SEXP sy,
                 SEXP sa, SEXP sb,
                 SEXP sX, SEXP sKmat,
                 SEXP sbeta0, SEXP stau20,
                 SEXP sknots, SEXP sdegree, SEXP sdiff,
                 SEXP sburnin, SEXP sstep,  SEXP siter)
{
  
  int burnin = Rf_asInteger(sburnin);
  int step   = Rf_asInteger(sstep);
  int iter   = Rf_asInteger(siter);
  int iterations = burnin + step*iter;
  
  int knots = Rf_asInteger(sknots);
  int degree = Rf_asInteger(sdegree);
  int diff = Rf_asInteger(sdiff);
  int n = Rf_nrows(sX);
  int K = Rf_ncols(sX); // knots + degree
  double *pX = REAL(sX);
  
  double tau2c = Rf_asReal(stau20);
  double aprime = Rf_asReal(sa) + (double)(K - diff)/2.0;
  double b = Rf_asReal(sb),  bprime;
  
  SEXP ssamp = PROTECT(allocMatrix(REALSXP, iter, K + 1 ));
  double *samp = REAL(ssamp);
  int acceptance = 0;
  
  double *betac = (double *) Calloc( K, double );
  double *betap = (double *) Calloc( K, double );
  
  double *beta0 = REAL(sbeta0);
  double *y = REAL(sy), *X = REAL(sX) ;
  
  /*variables de trabajo  */
  double *eta = (double *) Calloc( n, double );
  double mu;
  double *ytilde = (double *) Calloc( n, double );
  double *w = (double *) Calloc( n, double );
  double *mc = (double *) Calloc( K, double );
  double *mp = (double *) Calloc( K, double );
  double *temp = (double *) Calloc( K, double );
  
  SEXP Pc = PROTECT(allocMatrix(REALSXP, K , K ));
  SEXP Pp = PROTECT(allocMatrix(REALSXP, K , K ));
  
  double logold, lognew, qold, qnew, logdetPc, logdetPp, alpha;
  
  
  int i, j, index;
  
  Memcpy(betac, beta0, (size_t)K );
  
  //BLAS
  char *trans = "N";
  double blas_alpha=1.0, blas_beta=0.0;
  int inc=1;
  
  // Elementos de salida
  SEXP out = PROTECT(allocVector( VECSXP, 2));
  SEXP nm  = PROTECT(allocVector( STRSXP, 2));
  
  SET_STRING_ELT(nm, 0, mkChar("samp"));
  SET_STRING_ELT(nm, 1, mkChar("acceptance"));
  
  setAttrib(out, R_NamesSymbol, nm);
  
  Rprintf("Exp\naprime:%f\n", aprime);
  
  GetRNGstate();
  
  for( i =1; i <= iterations; i++ ){
    if( i > burnin && (i-burnin)%step == 0){
      index = (i-burnin)/step - 1;
      
      for ( j=0; j<K; j++)  samp[ iter*j + index ] = betac[j];
      
      samp[ iter*K + index ] = tau2c;
    }
    
    /* eta <- X betac */
    F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, pX, &n, 
             betac, &inc, &blas_beta, eta, &inc);
    
    logold = -0.5/tau2c*xSx( sKmat, betac );
    
    for(j=0; j < n ; j++ ){
      mu = logistic(eta[j]);
      logold += -log( mu ) - y[j]/mu;
      w[j] = (1.0 - mu)*(1.0 - mu);
      ytilde[j] = (y[j] - mu)/( mu*(1.0-mu) ) + eta[j];
    }
    
    tXWX( sX, w, 1.0/tau2c, sKmat, Pc);
    tXWytilde( sX, w, ytilde, mc ); 
    solve( Pc, mc);
    
    logdetPc = rmvnorm( mc, Pc, betap, "V"); //solo detSigma => "N"
    
    /* betap M-H step*/
    /* eta <- X betap */
    F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, pX, &n, 
             betap, &inc, &blas_beta, eta, &inc);
    
    lognew = -0.5/tau2c*xSx( sKmat, betap );
    
    for(j=0; j < n ; j++ ){
      mu = logistic(eta[j]);
      lognew += -log( mu ) - y[j]/mu;
      w[j] = (1.0 - mu)*(1.0 - mu);
      ytilde[j] = (y[j] - mu)/( mu*(1.0-mu) ) + eta[j];
    }
    
    tXWX( sX, w, 1.0/tau2c, sKmat, Pp);
    tXWytilde( sX, w, ytilde, mp ); 
    solve( Pp, mp);
    
    logdetPp = rmvnorm( NULL, Pp, NULL, "N"); //solo detSigma => "N"
    
    /*alpha*/
    for( j=0; j<K; j++) temp[j] = betap[j] - mc[j];
    qold = -0.5 *xSx( Pc, temp) + logdetPc;
    
    for( j=0; j<K; j++) temp[j] = betac[j] - mp[j];
    qnew = -0.5 *xSx( Pp, temp) + logdetPp;
    
    
    alpha = lognew + qnew - logold - qold;
    
    if( log(unif_rand()) < alpha  ){
      acceptance ++;
      Memcpy(betac, betap, (size_t)K );
      
    }
    
    /*Posterior samples for tau2c*/
    bprime = b + 0.5 * xSx( sKmat, betac );
    tau2c = 1.0/rgamma( aprime, 1.0/bprime );
    
    //Rprintf("ln:%f lo:%f qn:%f qo:%f xSx:%f\n", lognew, logold, qnew, qold, -0.5 *xSx( Pp, temp));
    //R_FlushConsole();
    R_CheckUserInterrupt();
  
  } // end i
  
  PutRNGstate();
  
  SET_VECTOR_ELT(out, 0, ssamp);
  SET_VECTOR_ELT(out, 1, ScalarReal( (double) acceptance / (double) iterations ));
  
  Free(betac);
  Free(betap);
  Free(eta);
  Free(ytilde);
  Free(w);
  Free(mc);
  Free(mp);
  Free(temp);
  
  UNPROTECT(5);
  return out;
}
