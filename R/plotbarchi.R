#' Plot barchi(t)
#'
#' @param fit List from tailindexreglogit function
#' @export

plotbarchi <- function(fit){
  link <- fit$link
  
  ylim <- if( link=="logit") c(-1,1) else c(-2,2) 
  link <- switch(link, "logit"=logistic, "exp"=exp )
  
  chibar <- 2*link(fit$etasamp) - 1
  mchibar <- apply(chibar, 1, mean)
  qchibar <- apply(chibar, 1, quantile,  probs =c(0.05, 0.95), names=FALSE )
  
  
  par(fig=c(0.57, 0.99, 0.25, 0.55), mar=c(0,3,0,0)+0.1)
	plot( fit$x, fit$y, pch=16, cex=0.5, axes=FALSE, col="gray",
	  xlab="", ylab="" )
	box()
	axis(2, at=1, cex=0.8)

	lines(fit$x, (mchibar+1)/2, col=1  ) 
	abline( h=mean(fit$y), col=2, lty=2)

  par(fig = c(0,1,0,1), mar=c(5,4,4,2)+0.1, new=TRUE )
  if( min(diff(fit$x)>=1) ){
    plot( as.Date(fit$x,origin="1970-01-01"), fit$y,  
          xlab="Date (in years)", ylab=expression(bar(chi)[t]), type="n", ylim=ylim )
  } else {
    plot( fit$x, fit$y, xlab="t", ylab=expression(bar(chi)[t]), ylim=ylim, type="n" )
  }

  polygon( c(fit$x,rev(fit$x)), c( qchibar[1,], rev(qchibar[2,])),
           density=NA, border=NA, col="gray") 
  lines( fit$x, mchibar , lty=1, lwd=2 )
  abline(h=c(-1,0,1), lty=2)
  abline(h=2*mean(fit$y)-1, col=2, lty=2)
  #if( !is.null(targetBarChi) )
  #  curve( targetBarChi(x), n=200, add=TRUE, col=2, lwd=2)
  



}

#' Plot chi(t)
#'
#' @param fit List from taildepreglogit function
#' @export
plotchi <- function(fit){

  chi <- logistic(fit$etasamp)
  mchi <- apply(chi, 1, mean)
  qchi <- apply(chi, 1, quantile,  probs =c(0.05, 0.95), names=FALSE )
   

  if( min(diff(fit$x)>=1) ){
    plot( as.Date(fit$x,origin="1970-01-01"), fit$y, pch="|", 
		xlab="Date (in years)", ylab=expression(chi[t]), type="p", ylim=c(0,1) )
  } else {
    plot( fit$x, fit$y, pch="|", xlab="t", ylab=expression(chi[t]), type="p", ylim=c(0,1))
  }

  polygon( c(fit$x,rev(fit$x)), c(qchi[1,], rev(qchi[2,])) , density=NA, border=NA, col="gray") 
  lines( fit$x, mchi , lty=1, lwd=2 )
  #if( !is.null(targetChi) )
  #  curve( targetChi(x), n=200, add=TRUE, col=2, lwd=2)


}
