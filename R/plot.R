#' Plot barchi
#' 
#' @param fit tailind output
#' @param xlab x-axis label
#' @param n Number of points for smooth line
#' 
#' @export
plotbarchi <- function( fit, xlab="t", ylim=c(-1,1), n=100){
  
  xx <- with(fit, seq(min(x), max(x), length.out = n))
  X <- with(fit, Bmat(xx, calcKnots(x, knots, degree), knots, degree ))
  
  h <- logistic
  if( fit$link == "log" ) h <- exp
  
  barchisamp <- 2*h(  X%*%t( fit$samp[,1:ncol(X)] )) - 1
  mbarchi <- apply(barchisamp, 1, mean)
  qbarchi <- apply(barchisamp, 1, quantile, probs=c(0.025, 0.975), names=FALSE)
  
  
  plot(xx, mbarchi, type="n", ylim=ylim,
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
  
  h <- logistic
  if( fit$link == "log" ) h <- exp
  
  etasamp <- h(  X%*%t( fit$samp[,1:ncol(X)] ))
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

#' Plot both measures chi_t and barchi_t
#' @export
plotMeasures <- function( fitdep, fitind, xlab="t", ylim=c(-1,1), n=100, target=NULL ){
  
  xx <- with(fitdep, seq(min(x), max(x), length.out = n))
  X <- with(fitdep, Bmat(xx, calcKnots(x, knots, degree), knots, degree ))
  
  chisamp <- logistic(  X%*%t( fitdep$samp[,1:ncol(X)] ))
  mchi <- apply(chisamp, 1, mean)
  qchi <- apply(chisamp, 1, quantile, probs=c(0.025, 0.975), names=FALSE)
  
  h <- logistic
  if( fitind$link == "log" ) h <- exp
  
  xx <- with(fitind, seq(min(x), max(x), length.out = n))
  X <- with(fitind, Bmat(xx, calcKnots(x, knots, degree), knots, degree ))
  
  barchisamp <- 2*h(  X%*%t( fitind$samp[,1:ncol(X)] )) - 1
  mbarchi <- apply(barchisamp, 1, mean)
  qbarchi <- apply(barchisamp, 1, quantile, probs=c(0.025, 0.975), names=FALSE)
  
  colp <- c("#c0c0c0", "#e0e0e0")
  coll <- c("red", "blue")
  ch <- rle(qbarchi[2,] > 1)
  
  if( length(ch$lengths) > 1 ){
    xpoly <- xx
    ypoly1 <- qchi[2,]
    ypoly2 <- qbarchi[2,]
    off <- 0
    for(i in 1:length(ch$lengths)){
      j <- ch$lengths[i]
      xpoly <- append(xpoly, rev(xx[off + (1:j)]), 2*off + j)
      ypoly1 <- append(ypoly1, rev(qchi[1,off + (1:j)]), 2*off + j)
      ypoly2 <- append(ypoly2, rev(qbarchi[1,off + (1:j)]), 2*off + j)
      off <- off + j
    }
    
    off <- 0
    for(i in 1:(length(ch$lengths)-1) ){
      j <- ch$lengths[i]
      xpoly <- append( xpoly, NA, off + 2*j ) 
      ypoly1 <- append( ypoly1, NA, off + 2*j )
      ypoly2 <- append( ypoly2, NA, off + 2*j )
      
      off <- off + 2*j + 1
    }  
  } else {
    xpoly <- c(xx, rev(xx))
    ypoly1 <- c( qchi[2,], rev(qchi[1,])  )
    ypoly2 <- c( qbarchi[2,], rev(qbarchi[1,])  )
  }
  
  par(mfrow=c(1,2))
  
  plot(fitdep$x, fitdep$y, pch="|",
       xlab=xlab, ylab=expression(chi[t]) )
  polygon(xpoly, ypoly1, col=colp[(!ch$values) + 1], border=NA )
  
  lines(xx, mchi)
  if( !is.null(target$dep) ){
    if( !is.matrix(target$dep) ) stop("target is not matrix")
    for( i in 2:ncol(target$dep) )
    lines(target$dep[,1], target$dep[,i], col=2, lty=i-1)
  }

  plot(xx, mbarchi, type="n", ylim=ylim,
       xlab=xlab, ylab=expression(bar(chi)[t]) )
  polygon(xpoly, ypoly2, col=colp[ch$values + 1], border=NA )
  
  lines(xx, mbarchi)
  abline( h=c(-1,0,1), lty=2 )
  
  
  
  if( !is.null(target$ind)){
    if( !is.matrix(target$ind) ) stop("target is not matrix")
    for( i in 2:ncol(target$ind) )
      lines(target$ind[,1], target$ind[,i], col=2, lty=i-1)
  }
  
}
