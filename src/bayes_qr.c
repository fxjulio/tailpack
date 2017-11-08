#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

double xSx( SEXP, double * );
double rgig( double, double, double );
void solve( SEXP A, double * );
double rmvnorm( double *, SEXP, double *, char *);

SEXP bqr(SEXP sbeta0, SEXP sX, SEXP sy, SEXP stau, SEXP stau20, SEXP sKdelta,
			SEXP sknots, SEXP sdegree, SEXP sdiff, 
			SEXP sburnin, SEXP sstep, SEXP siter,
			SEXP sa, SEXP sb, 
			SEXP sn0, SEXP ss0){

	int burnin = Rf_asInteger(sburnin);
	int step   = Rf_asInteger(sstep);
	int iter   = Rf_asInteger(siter);
	int iterations = burnin + step*iter;
	
	int knots = Rf_asInteger(sknots);
	int degree = Rf_asInteger(sdegree);
	int diff = Rf_asInteger(sdiff);
	
	int i, j, p, q, index;
	
	int k = Rf_ncols(sX);
	int n = Rf_nrows(sX);
	
	SEXP sbetac = PROTECT(allocVector(REALSXP, k));
	double *pbetac = REAL(sbetac);
	
	Memcpy(pbetac, REAL(sbeta0), (size_t)k );
	
	SEXP sPc = PROTECT(allocMatrix(REALSXP, k, k)); 
	SEXP smc = PROTECT(allocVector(REALSXP, k)  );
	double *pPc = REAL(sPc);
	double *pmc = REAL(smc);
	
	double tau2c = Rf_asReal(stau20);
	double atau2 = Rf_asReal(sa), btau2 =  Rf_asReal(sb);
	double aprime = atau2 + (double)(knots + degree - diff)/2.0;
	double bprime;
	
	double tau = Rf_asReal(stau);
	double mu = (1.0-2.0*tau)/(tau*(1.0-tau));
	double delta2 = 2.0/(tau*(1.0-tau));
	
	double sigmac, sast;
	double s0 = Rf_asReal(sn0);
	double nast = s0 + 3*n;
	
	double *vc = (double *) Calloc( n, double );
	double *temp = (double *) Calloc( n, double );
	
	double alpha, gamma;
	
	double *y = REAL(sy);
	double *pX = REAL(sX);
	double *pKdelta = REAL(sKdelta);
	
	//para BLAS
	char *trans = "N";
	double blas_alpha=1.0, blas_beta=0.0;
	int inc=1;
	
	SEXP ssamp = PROTECT(allocMatrix(REALSXP, iter, k+1 )); 
	double *psamp = REAL(ssamp); 
	
	//start iterations
	GetRNGstate();
	
	// init vc
	for( j=0; j<n; j++) vc[j] = exp_rand() ;
	
	
	Rprintf("bqr tau:%g, k:%d\n", tau, k);
	for( i=1; i<=iterations; i++){
		if( i > burnin && (i-burnin)%step == 0 ){
			index = (i-burnin)/step - 1;
			for( j=0; j < k; j++){
				psamp[iter*(j)+index] = pbetac[j];
			}
			psamp[iter*k + index] = tau2c;
		}
		
		/* Xbeta */
		F77_NAME(dgemv)(trans, &n, &k, &blas_alpha, pX, &n, 
			pbetac, &inc, &blas_beta, temp, &inc);
		
		/* posterior sigma */
		sast = s0;
		for( j=0; j<n; j++){
			sast += 2*vc[j] + R_pow_di( y[j] - temp[j] - mu*vc[j], 2 )/( delta2*vc[j] ) ;
		}		
		sigmac = 1.0/rgamma( nast/2.0, 1.0/(sast/2.0)  );
		
		/* posterior vc */
		gamma = 2.0/sigmac + mu*mu/(delta2*sigmac);
		for( j=0; j<n; j++){
			alpha = R_pow_di(y[j]-temp[j], 2)/( delta2*sigmac );
			vc[j] = rgig( 0.5, alpha, gamma );
		}
		
			
		/* Posterior betac */
		for( p=0; p < k; p++){
			for( q=p; q < k; q++){
				index = p*k + q;
				pPc[ index ]	= 1.0/tau2c * pKdelta[ index ];
				for( j=0; j < n; j++)
					pPc[ index ] +=  pX[ p*n + j ]*pX[ q*n + j ]/(delta2 * sigmac * vc[j]);
				
			}
		}
		for( p=0; p < k; p++){
			pmc[p] = 0.0;
			for( j=0; j < n; j++){
				pmc[ p ] +=  pX[ p*n + j ]*( y[j] - mu*vc[j] )/(delta2 * sigmac * vc[j]);
			}
		}
		
		//for( p=0; p < k; p++) pPc[ p*(k + 1) ] +=  1e-6;
		solve( sPc, pmc);
		rmvnorm( pmc, sPc, pbetac, "V");
		
		
		/*Posterior samples for tau2c*/
		bprime = btau2 + 1.0/2.0 * xSx( sKdelta, pbetac );
		tau2c = 1.0/rgamma( aprime, 1.0/bprime );
		
		R_CheckUserInterrupt();
		
	} // end i
	
	PutRNGstate();
	//end iterations
	

	Free(vc);
	Free(temp);
	UNPROTECT(4);
	return(ssamp);

}
