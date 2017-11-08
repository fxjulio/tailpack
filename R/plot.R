#' Plot barchi
#' 
#' @param fit tailind output
#' @param xlab x-axis label
#' @param n Number of points for smooth line
#' 
#' @export
plotbarchi <- function( fit, xlab="t", n=100){
  
  xx <- with(fit, seq(min(x), max(x), length.out = n))
  X <- with(fit, Bmat(xx, calcKnots(x, knots, degree), knots, degree ))
  
  barchisamp <- 2*logistic(  X%*%t( fit$samp[,1:ncol(X)] )) - 1
  mbarchi <- apply(barchisamp, 1, mean)
  qbarchi <- apply(barchisamp, 1, quantile, probs=c(0.025, 0.975), names=FALSE)
  
  
  plot(xx, mbarchi, type="n", ylim=c(-1,1),
       xlab=xlab, ylab=expression(bar(chi)[t]) )
  polygon(c(xx,rev(xx)), c(qbarchi[1,], rev(qbarchi[2,])), col="gray", border = NA)
  lines(xx, mbarchi)
  abline( h=c(-1,0,1), lty=2 )

}

#' Plot exponential samples
#' 
#' @param fit tailind output
#' @param xlab x-axis label
#' @param n Number of points for smooth line
#' 
#' @export
ploteta <- function( fit, xlab="t", n=100 ){
  
  xx <- with(fit, seq(min(x), max(x), length.out = n))
  X <- with(fit, Bmat(xx, calcKnots(x, knots, degree), knots, degree ))
  
  etasamp <- logistic(  X%*%t( fit$samp[,1:ncol(X)] ))
  meta <- apply(etasamp, 1, mean)
  qeta <- apply(etasamp, 1, quantile, probs=c(0.025, 0.975), names=FALSE)
  
  
  plot(fit$x, fit$y, pch=16, cex=0.6, col="gray", axes = FALSE,
       xlab=xlab, ylab=expression(eta[t]) )
  polygon(c(xx,rev(xx)), c(qeta[1,], rev(qeta[2,])), col="gray", border = NA)
  lines(xx, meta)
  axis(1)
  axis(2, at=1)
  box()
  
}

#' Plot chi
#' 
#' Plots a taildep ouput.
#' 
#' @param fit taildep output
#' @param xlab x-axis label
#' @param n Number of points for smooth line
#' @export
plotchi <- function( fit, xlab="t", n=100 ){
  
  xx <- with(fit, seq(min(x), max(x), length.out = n))
  X <- with(fit, Bmat(xx, calcKnots(x, knots, degree), knots, degree ))
  
  chisamp <- logistic(  X%*%t( fit$samp[,1:ncol(X)] ))
  mchi <- apply(chisamp, 1, mean)
  qchi <- apply(chisamp, 1, quantile, probs=c(0.025, 0.975), names=FALSE)
  
  
  plot(fit$x, fit$y, pch="|",
       xlab=xlab, ylab=expression(chi[t]) )
  polygon(c(xx,rev(xx)), c(qchi[1,], rev(qchi[2,])), col="gray", border = NA)
  lines(xx, mchi)
  
}
