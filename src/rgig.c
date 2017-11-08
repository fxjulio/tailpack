#include <R.h>
#include <Rmath.h>
#include <math.h>

double g( double, double, double );

double rgig( double lambda, double chi, double psi){
	double x, beta=sqrt(chi*psi), alpha=sqrt(chi/psi), m, x0, xast, k1, A1;
	double k2=0.0, A2=0.0, k3, A3, A, u, v, h;
	double xplus, lvplus, uplus;
	double a, b, c, p, q, phi, xminus, vplus, uminus;
	double aux1, aux2;
		
	//parameters lambda beta
	if( lambda > 2.0 || beta > 3.0  ){
		m = ( (lambda - 1.0) + sqrt((lambda - 1.0)*(lambda - 1.0) + beta*beta) ) / beta;
		a = -2.0*(lambda + 1.0)/beta - m;
		b = 2.0*(lambda - 1.0)/beta*m - 1.0;
		c = m;
		
		p = b - a*a/3.0;
		q = 2.0*a*a*a/27.0 - a*b/3.0 + c;
		
		phi = acos( - 0.5*q*sqrt(-27.0/(p*p*p)) );
		xminus = sqrt( -4.0/3.0*p )*cos( phi/3.0 + 4.0/3.0*M_PI ) -  a/3.0;
		xplus = sqrt( -4.0/3.0*p )*cos( phi/3.0 ) -  a/3.0;
		
		aux1 = 0.5*(lambda-1.);
		aux2 = 0.25*beta;
		
		lvplus = aux1*log(m) - aux2*(m+1./m);
		uminus = (xminus - m )*exp( aux1*log(xminus) - aux2*(xminus+1./xminus) - lvplus );
		uplus = (xplus - m)*exp( aux1*log(xplus) - aux2*(xplus+1./xplus) - lvplus ) ;
		
		do {
			u = uminus + unif_rand()*( uplus - uminus ) ;
			v = unif_rand();
			x = u/v + m;
			
		} while ( x<0 || log(v) > aux1*log(x) - aux2*(x + 1./x)  - lvplus );
		
		return x*alpha;
		
		
	}
	
	if ( lambda >= 1.0 - 2.25*beta*beta || beta > 0.2 ){
		
		m = beta/( (1.0-lambda) + sqrt( (1.0-lambda)*(1.0-lambda) + beta*beta ) );
		xplus = ( (1.0+lambda) + sqrt( (1.0+lambda)*(1.0+lambda) + beta*beta ) )/beta;
		
		vplus = sqrt( g(m, lambda, beta) );
		uplus = xplus * sqrt( g(xplus, lambda, beta) );
				
		do {
			u=unif_rand()*uplus;
			v=unif_rand()*vplus;
			x = u/v;
			
		} while ( v*v > g(x, lambda, beta) );
		
		return x*alpha; 
		
	}
	
	if( lambda >= 0. && beta > 0. ){
		
		// algorithm 1
		m = beta/( (1.0 - lambda) + sqrt((1.0-lambda)*(1.0-lambda) + beta*beta) );
		x0 = beta/( 1.0 - lambda );
		xast = fmax( x0, 2.0/beta );
		k1 = g(m, lambda, beta);
		A1 = k1*x0;
		if( x0 < 2.0/beta ){
			k2 = exp(-beta);
			if( lambda == 0.0 ){
				A2 = k2*log(2.0/(beta*beta));
			} else {
				A2 = k2*( pow(2.0/beta, lambda) - pow(x0, lambda) )/lambda; 
			}
		}
		k3 = pow(xast, lambda-1.0);
		A3 = 2*k3*exp( -xast*beta/2.0 )/beta;
		A = A1 + A2 + A3;
		
		do {
			u=unif_rand();
			v=unif_rand()*A;
			
			if( v <= A1 ){
				x = x0*v/A1;
				h = k1;
			} else if ( v <= A1 + A2) {
				v = v - A1;
				if( lambda == 0.0 )
					x = beta*exp( v*exp(beta) );
				else
					x = pow( pow(x0, lambda) + v*lambda/k2, 1.0/lambda );
				h = k2 * pow( x, lambda - 1.0 );
				
			} else {
				v = v - (A1 + A2);
				x = -2.0/beta*log( exp(-xast*beta/2.0) - v*beta/(2.0*k3) );
				h = k3*exp(-x*beta/2.0);
			}
			
		} while ( u*h > g(x, lambda, beta) );
		
		return x*alpha; 
		
	}
	
	
	return 0.0;
}

double g( double x, double lambda, double beta ){
	return pow(x, lambda - 1.0)*exp( -beta/2.0*( x + 1.0/x ) );
}
