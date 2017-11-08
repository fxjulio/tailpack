#' @title Tail index regression 2 (logit)
#'
#' @description
#' Posterior samples of barchi(t).
#' 
#' @param xy Matrix of data n rows and two columns.
#' @param t Time values. Vector of n.
#' @param threshold For quantile function. Default 0.95.
#' @param k Number of knots. Default 30.
#' @param l Degree of splines. Default 3.
#' @param d Diference between betas_k. Default 1.
#' @param atau2 Hyperparameter a for tau2. Default 0.001
#' @param btau2 Hyperparameter b for tau2. Default 0.001
#' @param step Thinning. Default 20
#' @param burnin First iterations discarted. Default 2000
#' @param iter Number of final iterations. Default 4000
#' @param lambda Penalization values. Default 1
#' @param maxiter Maximum of iterations for IWLS algorithm. Default 500
#'
#' @return A list
#' 
#' @export
tailindexreglogit2 <- function(xy, t, threshold=0.95, k=30, l=3, d=1,
  atau2=0.001, btau2=0.001, step=20, burnin=2000, iter=4000, 
  lambda=1, maxiter=500, control.bqr=list(burnin=1000, step=10, iter=1000), bqr=FALSE  ){
  
  if( class(t) == "Date" ){
    t <- as.numeric(t)
  }
  
  ############################
  # hacer variable exponencial
  
  z <- apply( xy, 1, min)
  
  if( !bqr ){
	u <- quantile( z, threshold,  names=FALSE )
	flag <- z > u
	y <- log( z[ flag ]/u  )
	bqrfit <- NA
  } else { #bqr
	  
	  bqrfit <- bqr( t, z, threshold, burnin = control.bqr$burnin, 
					step=control.bqr$step, iter=control.bqr$iter)
	  usamp <- bqrfit$X %*% t(bqrfit$samp[,1:ncol(bqrfit$X)])
	  u <- apply( usamp, 1, mean)
	  flag <- z > u
	  y <- log( z[flag]/u[flag] )
  }
    
  x <- t[flag]
  n <- sum(flag)


############################

  xi <- calcKnots( x, k, l )
  X <- Bmat(x, xi, k, l )
  Kmat <- Kmatrix( d, k, l)
  
###########################
# estimador IWLS
###########################
  stopifnot( mean(y) < 1 )
  b0 <- rep( logit(mean(y)), k+l)
  fit <- lapply( lambda, function(lambda) IWLSexpRint( X, Kmat, y, b0, lambda, maxiter, link="logit"))
  CVbar <- sapply( fit, function(x) ifelse(!x$convergence, x$CVbar, NA) )

  plot( lambda, CVbar, log="x" )

  fitopt <- fit[[which.min(CVbar)]]
  
########################
# valores iniciales
########################

beta0 <- fitopt$b
tau20 <- 10
  
  cat("lambdaopt:", fitopt$lambda, "\n")
  
########################
# Llamada a rutina C
########################
  iterations <- burnin + step*iter
  
  out <- .Call("bsmoothExp", y, atau2, btau2, X, Kmat, beta0, tau20,
				k, l, d, 
				burnin, step, iter)
    
    dimnames(out$samp)[[2]] <- c(paste0("beta",1:(k+l)), "tau2")

  list( x=x, y=y, u=u, z=z,
    samp = out$samp,
    acceptance = out$acceptance,
    CVbar=cbind(lambda=lambda, CVbar=CVbar), bqrfit=bqrfit
  )
}
