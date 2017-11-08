#' Bayesian quantile regression
#'
#' @export

bqr <- function(x, y, tau, knots=20, degree=3, diff=2,
			burnin=1000, step=10, iter=1000,
			n0=100, s0=1, a=0.001, b=0.001, tau20=10){
		
	xi <- calcKnots( x, knots, degree )
	X <- Bmat( x, xi, knots, degree )
	Kdelta <- Kmatrix(diff, knots, degree)
	
	b0 <- rep( quantile(y, tau, names=FALSE), ncol(X)  )
	
	samp <- .Call( "bqr",  b0=b0, X=X, y=y, tau=tau, tau20=tau20, Kdelta=Kdelta,
		knots=knots, degree=degree, diff=diff,
		burnin=burnin, step=step, iter=iter,
		a=a, b=b, n0=n0, s0=s0 )
	
	list( samp=samp, X=X, y=y, x=x )
}


