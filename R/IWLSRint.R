IWLSberRint<- function( X, Kmat, y, b0, lambda, maxiter=100, tol=1e-4, verbose=FALSE ){

	n <- length(y)
	
	out <- .Call("IWLS", X, Kmat, y, b0, lambda, maxiter, tol, as.numeric(verbose), "ber" )

	trH <- matrix.trace(out$H) ## se acerca a d(=1) cuando lambda es grande.

	list( b=out$b0, CVbar = sqrt(out$CV/n), trH=trH, lambda=lambda, convergence=out$convergence)
}

IWLSexpRint<- function( X, Kmat, y, b0, lambda, maxiter=100, tol=1e-4, verbose=FALSE, link="logit" ){

	n <- length(y)
  
  model <- paste0("exp", link)
  	
	out <- .Call("IWLS", X, Kmat, y, b0, lambda, maxiter, tol, as.numeric(verbose), model )

	trH <- matrix.trace(out$H) ## se acerca a d(=1) cuando lambda es grande.

	list( b=out$b0, nu=1.0, CVbar = sqrt(out$CV/n), trH=trH, lambda=lambda, convergence=out$convergence)
}



