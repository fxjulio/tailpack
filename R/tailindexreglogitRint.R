#' @title Tail index regression (logit)
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
#' @param Vtau2 Variance of tau2 samples. Default NULL
#' @param asigma2 Hyperparameter a for sigma2. Default 0.001
#' @param bsigma2 Hyperparameter b for sigma2. Default 0.001
#' @param anu Hyperparameter scale parameter. Default 0.1
#' @param bnu Hyperparameter scale parameter. Default 0.1
#' @param varnup Proposal parameter. Default 0.01
#' @param step Thinning. Default 20
#' @param burnin First iterations discarted. Default 2000
#' @param iter Number of final iterations. Default 4000
#' @param lambda Penalization values. Default 1
#' @param maxiter Maximum of iterations for IWLS algorithm. Default 500
#' @param link Character. Default "logit" or "exp" only.
#' @param theta0 Initial constant value. Default NULL.
#'
#' @return A list
#' 
#' @export
tailindexreglogit <- function(xy, t, threshold=0.95, k=30, l=3, d=1,
  atau2=0.001, btau2=0.001, Vtau2=NULL, asigma2=0.001, bsigma2=0.001,
  anu=0.1, bnu=0.1, varnup=0.01, step=10, burnin=2000, iter=4000, 
  lambda=1, maxiter=500, link="logit", theta0=NULL, control.bqr=NULL){
  
  if( class(t) == "Date" ){
    t <- as.numeric(t)
  }
  
  ############################
  # hacer variable exponencial
  
  z <- apply( xy, 1, min)
  
  if( is.null(control.bqr) ){
	u <- quantile( z, threshold,  names=FALSE )
	flag <- z > u
	y <- log( z[ flag ]/u  )
  } else { #bqr
	  
	  bqrfit <- bqr( t, z, threshold, burnin = 1000, step=10, iter=1000)
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
  
  U <- matrix(1, n, 1)

###########################
# estimador constante
###########################

  if( is.null(theta0) ){ #inicio gamma0
  
if( link == "logit" ){
loglik <- function(par, x, y){
	b0 <- par[1]
	nu <- par[2]
	eta <- b0
	mu <- logistic(eta)
	sum( dgamma(y, nu, scale=mu/nu, log=TRUE)   )
}

  if( mean(y) < 1 ){
    par <- c(logit(mean(y)), mean(y)^2/var(y))
    out <- optim( par, loglik, x=x, y=y, control=list(fnscale=-1) )
  } else {
    warning( "mean(y) > 1, theta0 high" )
    flagMean <- TRUE
  }

  
}

if( link == "exp" ){
loglik <- function(par, x, y){
	b0 <- par[1]
	nu <- par[2]
	eta <- b0
	mu <- exp(eta)
	sum( dgamma(y, nu, scale=mu/nu, log=TRUE)   )
}

  par <- c(log(mean(y)), mean(y)^2/var(y))
  out <- optim( par, loglik, x=x, y=y, control=list(fnscale=-1) )
}


if( link == "probit" ){
loglik <- function(par, x, y){
	b0 <- par[1]
	nu <- par[2]
	eta <- b0
	mu <- pnorm(eta)
	sum( dgamma(y, nu, scale=mu/nu, log=TRUE)   )
}

  par <- c(qnorm(mean(y)), mean(y)^2/var(y))
  out <- optim( par, loglik, x=x, y=y, control=list(fnscale=-1) )

}

###########################
# estimador IWLS
###########################

  b0 <- rep( out$par[1], k+l)
  nu0 <- out$par[2]
  fit <- lapply( lambda, function(lambda) IWLSexpRint( X, Kmat, y, b0, lambda, maxiter, link=link))
  CVbar <- sapply( fit, function(x) ifelse(!x$convergence, x$CVbar, NA) )

  plot( lambda, CVbar, log="x" )

  fitopt <- fit[[which.min(CVbar)]]
  
########################
# valores iniciales
########################

intercept <- bsIntercept(fitopt$b, k, l)

gamma0 <- intercept
beta0 <- fitopt$b - gamma0
nu0 <- fitopt$nu

} # end gamma0
else {
  gamma0 <- theta0
  beta0 <- rep(0, k+l)
  nu0 <- 1
  fitopt <- list(lambda=10)
}
  
  if( !is.null(Vtau2) ){
    Etau2 <- 1/fitopt$lambda
    atau2 <- 2 + Etau2^2 / Vtau2
    btau2 <- (atau2 - 1)*Etau2
    cat("atau2:", atau2, "btau2:", btau2, "\n")
  }
  
  cat("lambdaopt:", fitopt$lambda ,"gamma0:", gamma0, "nu0:", nu0, "\n")
  
########################
# Llamada a rutina C
########################
  iterations <- burnin + step*iter
  
  out <- .Call("tail", X, Kmat, k, l, d,
    beta0, gamma0, nu0, y,
    atau2, btau2, asigma2, bsigma2,
    anu, bnu, varnup,
    step, burnin, iter, "gamma", link)
    
    dimnames(out$betasamp)[[1]] <- paste0("beta",1:(k+l))

  list( x=x,
    y=y, u=u, z=z,
    etasamp=X%*%out$betasamp + U%*%matrix(out$gammasamp, 1, iter),
    samp = with(out, cbind(tau2=tau2samp, sigma2=sigma2samp, gamma=gammasamp, nu=nusamp, beta=t(betasamp))),
    acceptance=with(out, c(beta=acceptance_beta, gamma=acceptance_gamma, nu=acceptance_nu)/iterations),
    priors=c(atau2=atau2, btau2=btau2, asigma2=asigma2, bsigma2=bsigma2),
    CVbar=cbind(lambda=lambda, CVbar=CVbar), link=link, bqrfit=bqrfit
  )
}
