#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <math.h>

static union{
  double d;
  struct {
    int j, i;
  } n;
} eco;

#define EXP_A (1048576/M_LN2)
#define EXP_C 60801

#define EXP(y) ( eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), eco.d )



//R CMD config
//R CMD SHLIB -LC:/ARCHIV~1/R/R-30~1.1/bin/i386 -lRblas myfunctions.c
//R CMD SHLIB -LC:/DOCUME~1/jeavila/MISDOC~1/R/R-31~1.3/bin/i386 -lRlapack -LC:/ARCHIV~1/R/R-30~1.1/bin/i386 -lRblas myfunctions.c

//SEXP getListElement(SEXP, const char *);
double xSx( SEXP , double *  );
void tXWX( SEXP , double *, double, SEXP, SEXP  );
void tXWytilde( SEXP, double *, double *, double * );
double loglikBerLogit( SEXP, double *, double *, double *, double *, double );
double loglikExpLogit( SEXP, double *, double *, double *, double *, double );
double loglikExpExp( SEXP, double *, double *, double *, double *, double );
double loglikExpProbit( SEXP, double *, double *, double *, double *, double );
void solve( SEXP A, double * );
double rmvnorm( double *, SEXP, double *, char *);
double bsIntercept( double *, int, int );
double logliknu( SEXP , double *, double );
double logistic(double);

typedef double (*LOGLIK)(SEXP, double *, double *, double *, double *, double );
 // LOGLIK( y, *etatilde, *eta, *w, *ytilde, nu)
 
typedef double (*LINK)(double);

SEXP tail( SEXP X, SEXP Kmat, SEXP sk, SEXP sl, SEXP sd, 
  SEXP beta0, SEXP gamma0, SEXP nu0, SEXP y, 
  SEXP satau2, SEXP sbtau2,
  SEXP sasigma2, SEXP sbsigma2,
  SEXP sanu, SEXP sbnu, SEXP svarnup,
  SEXP sstep, SEXP sburnin,  SEXP siter, SEXP model, SEXP link ){

  double atau2 = Rf_asReal(satau2);
  double btau2 = Rf_asReal(sbtau2);
  double asigma2 = Rf_asReal(sasigma2);
  double bsigma2 = Rf_asReal(sbsigma2);
  double varnup = Rf_asReal(svarnup);
  int step   = Rf_asInteger(sstep);
  int burnin = Rf_asInteger(sburnin);
  int iter   = Rf_asInteger(siter);
  int k = Rf_asInteger(sk);
  int l = Rf_asInteger(sl);
  int d = Rf_asInteger(sd);
  int n = Rf_nrows(X);
  int K = Rf_ncols(X); // k+l
  SEXP betasamp = PROTECT(allocMatrix(REALSXP, K , iter ));
  SEXP tau2samp = PROTECT(allocVector(REALSXP, iter)  );
  SEXP gammasamp = PROTECT(allocVector(REALSXP, iter) );
  SEXP sigma2samp = PROTECT(allocVector(REALSXP, iter)  );
  SEXP Pc = PROTECT(allocMatrix(REALSXP, K , K ));
  SEXP Pp = PROTECT(allocMatrix(REALSXP, K , K ));
  double *rbetasamp = REAL( betasamp );
  double *ptau2samp = REAL( tau2samp );
  double *psigma2samp = REAL( sigma2samp );
  double *pgammasamp = REAL( gammasamp );
  double *rX = REAL(X);
  
  int acceptance_beta = 0;
  int acceptance_gamma = 0;
  int acceptance_nu = 0;
  
  double *betac = (double *) Calloc( K, double );
  double *betap = (double *) Calloc( K, double );
  double tau2c = 10.0, sigma2c = 10.0;
  double *etac = (double *) Calloc( n, double );
  double *etatilde = (double *) Calloc( n, double );
  double *ytilde = (double *) Calloc( n, double );
  double *w = (double *) Calloc( n, double );
  double *mc = (double *) Calloc( K, double );
  double *betamode = (double *) Calloc( K, double );
  double *mp = (double *) Calloc( K, double );
  double *temp = (double *) Calloc( K, double );
  double intercept, gammac, gammap;
  double nuc = Rf_asReal(nu0), nup;
  
  int iterations = burnin + step * iter;
  int index, i, j;
  // int jj;
  
  // Elementos de salida
  SEXP out = PROTECT(allocVector( VECSXP, 8));
  SEXP nm  = PROTECT(allocVector( STRSXP, 8));
  
  SET_STRING_ELT(nm, 0, mkChar("betasamp"));
  SET_STRING_ELT(nm, 1, mkChar("tau2samp"));
  SET_STRING_ELT(nm, 2, mkChar("gammasamp"));
  SET_STRING_ELT(nm, 3, mkChar("sigma2samp"));
  SET_STRING_ELT(nm, 4, mkChar("nusamp"));
  SET_STRING_ELT(nm, 5, mkChar("acceptance_beta"));
  SET_STRING_ELT(nm, 6, mkChar("acceptance_gamma"));
  SET_STRING_ELT(nm, 7, mkChar("acceptance_nu"));
  setAttrib(out, R_NamesSymbol, nm);
  
  double logold, lognew;
  double logdetPc, logdetPp;
  double qold, qnew, alpha;
  double aprime, bprime;
  double gamma_max = 6.0;
  
  LOGLIK loglik;
  int modExp = 0;
  double anup, scalenup;
  double anu = Rf_asReal(sanu), bnu = Rf_asReal(sbnu);
  
  SEXP nusamp = PROTECT(allocVector(REALSXP, iter) );
  double *pnusamp = REAL( nusamp );
  
  if( strcmp( CHAR(asChar(model)), "ber") == 0  && 
      strcmp( CHAR(asChar(link)), "logit") == 0 ){
  
    loglik = loglikBerLogit;
    modExp = 0;
    Rprintf("ber logit\n");  
  } else if( strcmp( CHAR(asChar(model)), "gamma") == 0 && 
      strcmp( CHAR(asChar(link)), "logit") == 0   ){
    loglik = loglikExpLogit;
    modExp = 1;
    Rprintf("gamma logit\n");  
  } else if( strcmp( CHAR(asChar(model)), "gamma") == 0 && 
      strcmp( CHAR(asChar(link)), "exp") == 0   ){
    loglik = loglikExpExp;
    modExp = 1;
    Rprintf("gamma exp\n");  
  } else if( strcmp( CHAR(asChar(model)), "gamma") == 0 && 
      strcmp( CHAR(asChar(link)), "probit") == 0   ){
    loglik = loglikExpProbit;
    modExp = 1;
    Rprintf("gamma probit\n");  
  } else {
    error("model or link not adequate\n");
  }
  
  
  R_FlushConsole();
  
  //para BLAS
  char *trans = "N";
  double blas_alpha=1.0, blas_beta=0.0;
  int inc=1;
  
  // valores iniciales
  gammac = Rf_asReal(gamma0);
  F77_NAME(dcopy)(&K, REAL(beta0), &inc, betac, &inc); // beta0 has zero mean
  F77_NAME(dcopy)(&K, REAL(beta0), &inc, betamode, &inc); // betamode
  
  for( j=0; j<n; j++) etac[j] = gammac;
  blas_beta = 1.0;
  F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
    betac, &inc, &blas_beta, etac, &inc);
    
  GetRNGstate();
  
  gammap = 0;
  //Rprintf("logdetPc:%f logdetPp:%f", logdetPc, logdetPp);
  
  
  for( i =1; i <= iterations; i++ ){
    if( i > burnin && (i-burnin)%step == 0){
      index = (i-burnin)/step - 1;
      F77_NAME(dcopy)( &K, betac, &inc, rbetasamp+K*index, &inc );
      ptau2samp[index] = tau2c;
      psigma2samp[index] = sigma2c;
      pgammasamp[index] = gammac;
      if( modExp ) pnusamp[index] = nuc;
      //Rprintf("\b\b\b\b\b%03.1f", (double) i/iterations*100);
      //R_FlushConsole();
    }
    
    
    // for( jj=1; jj <= 15 && 1 <= i && i <= 200; jj++){
      // blas_beta = 0.0;
      // F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
        // betac, &inc, &blas_beta, etatilde, &inc);
      // loglik( y, etatilde, etac, w, ytilde, nuc );
      // tXWX( X, w, 1.0/tau2c, Kmat, Pc); // pasos WLS
      // tXWytilde( X, w, ytilde, mc );
      // solve( Pc, mc);
      
      // for( j=0; j<K; j++) temp[j] = mc[j] - betac[j];
      // blas_beta = 1.0;
      // F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
        // temp, &inc, &blas_beta, etac, &inc);
        
      // intercept = bsIntercept( mc, k, l);
      // for( j=0; j<K; j++) mc[j] -= intercept;
      // gammac += intercept;
        
      // F77_NAME(dcopy)( &K, mc, &inc, betac, &inc );
      
      // //Rprintf("jj:%d\n", jj);
      
    // }
    
    
    /* posterior samples beta */
    /* usando betac */
    
    blas_beta = 0.0;
    F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
      betac, &inc, &blas_beta, etatilde, &inc); // etatilde: Xbetac
    logold = -0.5 *xSx( Kmat, betac)/tau2c;
    logold += loglik( y, etatilde, etac, w, ytilde, nuc );
    
    /* intercambiando betac por betamode en etac*/
    // for( j=0; j<K; j++) temp[j] = betamode[j] - betac[j];
    // blas_beta = 1.0;
    // F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
      // temp, &inc, &blas_beta, etac, &inc);
    // /*calculando etatilde = X*betamode  */
    // blas_beta = 0.0;
    // F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
      // betamode, &inc, &blas_beta, etatilde, &inc); // etatilde = X*betamode
    
    /*calculando w , ytilde en base a betamode */    
    loglik( y, etatilde, etac, w, ytilde, nuc );
    
    tXWX( X, w, 1.0/tau2c, Kmat, Pc);
    
    tXWytilde( X, w, ytilde, mc ); // aqui betamode se sobreescribe mc <- betamode
    
    solve( Pc, mc);

    logdetPc = rmvnorm( mc, Pc, betap, "V"); //solo detSigma => "N"
    //rmvnorm( mc, Pc, betap, "V");

    /*intercambiando betamode por betap*/
    //for( j=0; j<K; j++) temp[j] = betap[j] - betamode[j];
    for( j=0; j<K; j++) temp[j] = betap[j] - betac[j];
    blas_beta = 1.0;
    F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
      temp, &inc, &blas_beta, etac, &inc);
    
    blas_beta = 0.0;
    F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
     betap, &inc, &blas_beta, etatilde, &inc);     // esto no es necesario // separar loglik de do_W_ytilde
     
    intercept = bsIntercept( betap, k, l);
    lognew = gammac + intercept <= 6 ? -0.5 *xSx( Kmat, betap)/tau2c : -DBL_MAX;
    lognew += loglik( y, etatilde, etac, w, ytilde, nuc);
    
    tXWX( X, w, 1.0/tau2c, Kmat, Pp); // esto no se usará
    
    tXWytilde( X, w, ytilde, mp );// esto no se usará
    solve(Pp, mp);// esto no se usará
        
    logdetPp = rmvnorm( NULL, Pp, NULL, "N"); //solo logdetP => "N"// esto no se usará
    
    for( j=0; j<K; j++) temp[j] = betap[j] - mc[j];
    // qold = -0.5 *xSx( Pc, temp); // forma betamode
    
    qold = -0.5 *xSx( Pc, temp) + logdetPc;  //forma antigua
    
    for( j=0; j<K; j++) temp[j] = betac[j] - mp[j]; //forma antigua
    qnew = -0.5 *xSx( Pp, temp) + logdetPp;
    
    // for( j=0; j<K; j++) temp[j] = betac[j] - mc[j];
    // qnew = -0.5 *xSx( Pc, temp);
    
    /*dejar betamode <- mc */
    // intercept = bsIntercept( mc, k, l);
    // F77_NAME(dcopy)( &K, mc, &inc, betamode, &inc );
    // for( j=0; j<K; j++) betamode[j] -= intercept;
    
    
    alpha = lognew + qnew - logold - qold;
    
    //Rprintf("\b\b\b\b\b\bgc:%f int:%f \n", gammac, bsIntercept( betap, k, l) );
    
    if( log(unif_rand()) < alpha  ){
      acceptance_beta ++;
      intercept = bsIntercept( betap, k, l);
      gammac += intercept;
      for( j=0; j<K; j++) betap[j] -= intercept;
      F77_NAME(dcopy)( &K, betap, &inc, betac, &inc );
      
    } else {
      for( j=0; j<K; j++) temp[j] = betac[j] - betap[j];
      blas_beta = 1.0;
      F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, rX, &n, 
        temp, &inc, &blas_beta, etac, &inc);
    
    }
    
    /*posterior samples for gammac*/
    for( j=0; j < n; j++) etatilde[j] = gammac;
    logold = gammac <= gamma_max ? -0.5*gammac*gammac/sigma2c: -DBL_MAX;
    logold += loglik( y, etatilde, etac, w, ytilde, nuc );
    REAL(Pc)[0] = 0.0; // solo usar la posicion 1
    mc[0] = 0.0;
    for( j=0; j < n; j++){
      REAL(Pc)[0] += w[j];
      mc[0] += w[j]*ytilde[j];
    }
    REAL(Pc)[0] += 1.0/sigma2c; 
    mc[0] = mc[0]/REAL(Pc)[0];
    
    gammap = mc[0] + 1.0/sqrt(REAL(Pc)[0])*rnorm(0,1);
    
    for( j=0; j < n; j++){
      etatilde[j] = gammap;
      etac[j] += gammap - gammac;
    }
    lognew = gammap <= gamma_max ? -0.5*gammap*gammap/sigma2c :  -DBL_MAX;
    lognew += loglik( y, etatilde, etac, w, ytilde, nuc );
    REAL(Pp)[0] = 0.0; // solo usar la posicion 1
    mp[0] = 0.0;
    for( j=0; j < n; j++){
      REAL(Pp)[0] += w[j];
      mp[0] += w[j]*ytilde[j];
    }
    REAL(Pp)[0] += 1.0/sigma2c; 
    mp[0] = mp[0]/REAL(Pp)[0];
    
    qold = -0.5 * REAL(Pc)[0]*(gammap - mc[0])*(gammap - mc[0]) + log(REAL(Pc)[0]);
    qnew = -0.5 * REAL(Pp)[0]*(gammac - mp[0])*(gammac - mp[0]) + log(REAL(Pp)[0]);

    alpha = lognew + qnew - logold - qold;
    
    //Rprintf("gc:%f gp:%f a:%e\n", gammac, gammap, alpha);
    
    if( log(unif_rand()) < alpha  ){
      acceptance_gamma ++ ;
      gammac = gammap;
    } else {
      for( j=0; j < n; j++) etac[j] += gammac - gammap;
    }
    
    
    /*Posterior samples for nuc */
    if( modExp && 0 ){
      anup = nuc*nuc / varnup;
      scalenup = varnup/nuc; 
      logold = logliknu( y, etac, nuc ) + dgamma( nuc, anu, 1.0/bnu, 1);
      nup = rgamma( anup, scalenup);
      lognew = logliknu( y, etac, nup ) + dgamma( nup, anu, 1.0/bnu, 1);
      
      qold = dgamma( nup, anup, scalenup, 1);
      
      anup = nup*nup / varnup;
      scalenup = varnup/nup; 
      qnew = dgamma( nuc, anup, scalenup, 1);
            
      alpha = lognew + qnew - logold - qold;
    
      if( log(unif_rand()) < alpha  ){
        nuc = nup;
        acceptance_nu++;
      }
    
    }
    
    
    /*Posterior samples for tau2c*/
    aprime = atau2 + (double)(K - d)/2.0;
    bprime = btau2 + 1.0/2.0 * xSx( Kmat, betac );
    
    tau2c = 1.0/rgamma( aprime, 1.0/bprime );
    
    /*Posterior samples for sigma2c*/
    aprime = asigma2 + 0.5;
    bprime = bsigma2 + 1.0/2.0 * gammac*gammac;
    
    sigma2c = 1.0/rgamma( aprime, 1.0/bprime );
    
    R_CheckUserInterrupt();
        
  }
  Rprintf("\n");
  
  PutRNGstate();
  
  Free(betac);
  Free(betap);
  Free(etac);
  Free(ytilde);
  Free(w);
  Free(mc);
  Free(mp);
  Free(temp);
  
  SET_VECTOR_ELT(out, 0, betasamp);
  SET_VECTOR_ELT(out, 1, tau2samp);
  SET_VECTOR_ELT(out, 2, gammasamp);
  SET_VECTOR_ELT(out, 3, sigma2samp);
  SET_VECTOR_ELT(out, 4, nusamp);
  SET_VECTOR_ELT(out, 5, ScalarReal((double)acceptance_beta));
  SET_VECTOR_ELT(out, 6, ScalarReal((double)acceptance_gamma));
  SET_VECTOR_ELT(out, 7, ScalarReal((double)acceptance_nu));
  
  UNPROTECT(9);
  return out;

 }
 
void solveMat( SEXP A, SEXP B ){
  int n = Rf_nrows(A);
  int *ipiv = (int *) Calloc(n, int);
  int K = Rf_ncols(B);
  char *uplo = "L";
  
  double *avals = (double *) Calloc( n * n, double);
  Memcpy(avals, REAL(A), (size_t)n * n);
  int info, lwork = n+1;
  double *work = (double *) Calloc( lwork, double);
 
  F77_CALL(dsysv)(uplo, &n, &K, avals, &n, ipiv, REAL(B), &n, work, &lwork, &info);
  
  Free(ipiv);
  Free(avals);
  Free(work);
  
  if( info < 0 )
    error("argument %d of Lapack routine %s had invalid value",
        -info, "dgesv");
  if (info > 0)  
    error("Lapack routine %s: system is exactly singular: U[%d,%d] = 0",
      "dgesv", info, info);
      
}


void logistic2(double *x, double *ans ){
  //double ans = 1.0;
  double a[5];
  int i;
  a[0] = 1.0/12.0;
  a[1] = 1.0/10.0;
  a[2] = 17.0/168.0;
  a[3] = 31.0/306.0;
  a[4] = 0.1013196481;
  
  *ans = 1.0;
  for (i = 4; i >= 0; --i )
    *ans = 1.0 - (*x)*(*x) * (*ans) * a[i];
    
  *ans = 1.0/2.0 + (*x)*(*ans)/4.0 ;

}


// etatilde es Xbeta o Ugamma
double loglikBerLogit( SEXP y, double *etatilde, double *eta, double *w, double *ytilde, double nu ){
  double ll = 0.0;
  double *py = REAL(y), mu, v;
  int n = LENGTH(y);
  
  for( int i=0; i < n; i++ ){
      mu = logistic(eta[i]);
      ll += py[i]*log(mu) + (1.0-py[i])*log(1.0-mu);
      v = mu * (1.0 - mu);
      ytilde[i] = (py[i] - mu)/v + etatilde[i];
      w[i] = v;
  }

  return ll;
}


// etatilde es Xbeta o Ugamma
double loglikExpLogit( SEXP y, double *etatilde, double *eta, double *w, double *ytilde, double nu){
  double ll = 0.0;
  double *py = REAL(y), mu;
  int n = LENGTH(y);
  
  for( int i=0; i < n; i++ ){
      mu = logistic(eta[i]);
      ll += -nu*log(mu) + ( nu-1.0 )*log(py[i]) - nu * py[i]/mu;
      ytilde[i] = (py[i] - mu)/(mu * (1.0 - mu)) + etatilde[i];  
      w[i] = nu*(1-mu)*(1-mu);
  }

  return ll;
}

double loglikExpExp( SEXP y, double *etatilde, double *eta, double *w, double *ytilde, double nu){
  double ll = 0.0;
  double *py = REAL(y), mu;
  int n = LENGTH(y);
  
  for( int i=0; i < n; i++ ){
      mu = exp(eta[i]);
      ll += -nu*log(mu) + ( nu-1.0 )*log(py[i]) - nu * py[i]/mu;
      ytilde[i] = (py[i] - mu)/mu + etatilde[i];  
      w[i] = nu;
  }

  return ll;
}

// etatilde es Xbeta o Ugamma
double loglikExpProbit( SEXP y, double *etatilde, double *eta, double *w, double *ytilde, double nu){
  double ll = 0.0;
  double *py = REAL(y), mu, etatemp, gprime;
  int n = LENGTH(y);
  
  for( int i=0; i < n; i++ ){
      etatemp =  (-6.0 <= eta[i] && eta[i] <= 6.0) ? eta[i] : 6.0;
      if( eta[i] < -6.0 ) etatemp = -6.0;
      mu = pnorm(etatemp, 0, 1, 1, 0);
      ll += -nu*log(mu) + ( nu-1.0 )*log(py[i]) - nu * py[i]/mu;
      gprime = dnorm(etatemp, 0, 1, 0);
      ytilde[i] = (py[i] - mu)/gprime + etatilde[i];  
      w[i] = nu*gprime*gprime/(mu*mu);
  }

  return ll;
}


double logliknu( SEXP y, double *eta, double nu){
  double ll = 0.0;
  double *py = REAL(y), mu;
  int n = LENGTH(y);
  
  for( int i=0; i < n; i++ ){
      mu = logistic(eta[i]);
      ll += dgamma( py[i], nu, mu/nu, 1 );
      //ll += nu*log(nu) - lgamma(nu) -nu*log(mu) + ( nu-1.0 )*log(py[i]) - nu * py[i]/mu;
  }

  return ll;
}


double bsIntercept( double *b, int k, int l ){
  int p = k + l;
  double *w = (double *)Calloc( p, double );
  if( l == 2 ){
    for( int i=0; i < p ; i++ ) w[i] = 1.0;
    w[0] = 1.0/6.0;
    w[p-1] = 1.0/6.0;
    w[1] = 5.0/6.0;
    w[p-2] = 5.0/6.0;
  }
  
  if( l == 3 ){
    for( int i=0; i < p ; i++ ) w[i] = 1.0;
    w[0] = 1.0/24.0;
    w[p-1] = 1.0/24.0;
    w[1] = 12.0/24.0;
    w[p-2] = 12.0/24.0;
    w[2] = 23.0/24.0;
    w[p-3] = 23.0/24.0;
  }
  
  double sum = 0.0;
  for( int i=0; i < p ; i++ ) sum += b[i] * w[i];
  
  Free(w);
  
  return sum / (double) k ;

}

SEXP IWLS(SEXP sX, SEXP sKmat, SEXP sy, SEXP sb0, SEXP slambda,
  SEXP smaxiter, SEXP stol, SEXP sverbose, SEXP smodel){

  double lambda = Rf_asReal(slambda);
  double tol = Rf_asReal(stol);
  int maxiter = Rf_asInteger(smaxiter);
  int verbose = Rf_asInteger(sverbose);
  int convergence = 0;
  double crit = 1, dtemp;
  int iter = 1, i, j;
  int n = Rf_nrows(sX);
  int K = Rf_ncols(sX);
  
  SEXP b0 = PROTECT(allocVector(REALSXP, K ));
  SEXP b1 = PROTECT(allocVector(REALSXP, K ));
  SEXP P0 = PROTECT(allocMatrix(REALSXP, K, K) );
  SEXP H = PROTECT(allocMatrix(REALSXP, n, n) );
  SEXP mTemp = PROTECT(allocMatrix(REALSXP, K, n) );
  
  double *pb0 = REAL(b0);
  double *pb1 = REAL(b1);
  double *pTemp = REAL(mTemp);
  double *pX = REAL(sX);
  double *pb00 = REAL(sb0);
  
  LOGLIK loglik;
  LINK link;
  
  if( strcmp( CHAR(asChar(smodel)), "ber") == 0 ){
    loglik = loglikBerLogit;
    link = logistic;
  } else if ( strcmp( CHAR(asChar(smodel)), "explogit") == 0  ){
    loglik = loglikExpLogit;
    link = logistic;
  } else if ( strcmp( CHAR(asChar(smodel)), "explog") == 0  ){
    loglik = loglikExpExp;
    link = exp;
  } else {
    error("Model is not \"ber\" or \"explogit\" \"expexp\"\n");
  }
  
  for( i=0; i < K; i++) pb0[i] = pb00[i];
  
  // Elementos de salida
  SEXP out = PROTECT(allocVector( VECSXP, 4));
  SEXP nm  = PROTECT(allocVector( STRSXP, 4));
  
  SET_STRING_ELT(nm, 0, mkChar("b0"));
  SET_STRING_ELT(nm, 1, mkChar("CV"));
  SET_STRING_ELT(nm, 2, mkChar("H"));
  SET_STRING_ELT(nm, 3, mkChar("convergence"));
  setAttrib(out, R_NamesSymbol, nm);
  
  //para BLAS
  char *trans = "N";
  double blas_alpha=1.0, blas_beta=0.0;
  int inc=1;
  
  double *eta0 = (double *) Calloc( n, double );
  double *w0 = (double *) Calloc( n, double );
  double *ytilde = (double *) Calloc( n, double );
  
  while ( crit > tol ){

    F77_NAME(dgemv)(trans, &n, &K, &blas_alpha, pX, &n, 
      pb0, &inc, &blas_beta, eta0, &inc);

    loglik( sy, eta0, eta0, w0, ytilde, 1.0 );
    
    tXWX( sX, w0, lambda, sKmat, P0);
    
    tXWytilde( sX, w0, ytilde, pb1 ); 
    
    solve( P0, pb1);
      
    crit = 0.0;
    dtemp = 0.0;
    for( i=0; i < K ; i++){
      crit += R_pow_di( pb1[i]-pb0[i], 2);
      dtemp += R_pow_di( pb0[i], 2);
      pb0[i] = pb1[i];
    }
    crit /= dtemp;
    
    if( verbose ) Rprintf("iter:%d crit:%e\n", iter, crit);

    iter ++;
    if( iter > maxiter ){
      convergence = 1;
      break;
    }  
  }
  
  // H = X (X'WX + lambda Kmat)^-1 X'W
  
  for( i=0; i < n; i++ ){
    for( j=0; j < K; j++)   pTemp[ K*i + j ] = pX[ n*j + i ] * w0[i];
  }
  
  solveMat( P0, mTemp  );
  
  
  blas_beta = 0.0;
  F77_NAME(dgemm)(trans, trans, &n, &n, &K, &blas_alpha, pX, &n, REAL(mTemp), &K, &blas_beta, REAL(H), &n);
  
  
  //CV
  double CV = 0.0, mu0;
  for( i=0; i < n; i++){
    mu0 = link(eta0[i]);
    CV += R_pow_di(REAL(sy)[i]-mu0, 2) / R_pow_di(1.0 - REAL(H)[(n+1)*i], 2);
  }
  
  SET_VECTOR_ELT(out, 0, b0);
  SET_VECTOR_ELT(out, 1, ScalarReal(CV));
  SET_VECTOR_ELT(out, 2, H);
  SET_VECTOR_ELT(out, 3, ScalarReal((double)convergence));
  
  Free(eta0);
  Free(w0);
  Free(ytilde);
  UNPROTECT(7);
  
  return out;
  
}



// /* get the list element named str, or return NULL */

// SEXP getListElement(SEXP list, const char *str){
    // SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

    // for (int i = 0; i < length(list); i++)
      // if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
         // elmt = VECTOR_ELT(list, i);
         // break;
      // }
    // return elmt;
// }

SEXP myexp( SEXP sx ){
  double x = Rf_asReal(sx);
  
  return ScalarReal( (double)EXP(x) ) ;

}

