#' Plot copula chi(u)
#' @examples 
#' blevd <- function(u,v, ... ) exp(-( (-log(u))^(1/alpha)+(-log(v))^(1/alpha))^alpha )
#' alpha <- c(0.1,0.3,0.5,0.7,0.9)
#' plotCopulaChi(blevd, u=seq(0.9,0.9999, length.out = 100),
#'               alpha=alpha  )
#' points(rep(1, length(alpha)), 2-2^alpha, col=2)
#' @export
plotCopulaChi <- function( copula, u=seq(0.9,0.9999, length.out = 100), ...){
  
  chiu <- mapply(function( u, ... ){
    2 - log( copula( u, u, ... ) )/log( u ) # sec 3.3.1 Coles1999
  }, u=u, MoreArgs = list(...) )
  

  if( is.vector(chiu)){
    plot(u, chiu, xlab="u", ylab=expression(chi(u)))
  }  else {
    plot(0,0, xlab="u", ylab=expression(chi(u)),
         type = "n", xlim = range(u), ylim = range(chiu))
    for(i in 1:nrow(chiu)) lines(u, chiu[i,])
  }

}


#' Plot copula barchi
#' @examples
#' blevd <- function(u,v, ... ) exp(-( (-log(u))^(1/alpha)+(-log(v))^(1/alpha))^alpha )
#' alpha <- c(0.1,0.3,0.5,0.7,0.9)
#' plotCopulaBarChi(blevd, u=seq(0.9,0.9999, length.out = 100),
#'                  alpha=alpha  )
#' points(1, 1, col=2)
#' @export
plotCopulaBarChi <- function( copula, u=seq(0.9,0.9999, length.out = 100), ...){
  
  barchiu <- mapply(function( u, ... ){
    2*log(1-u)/log( 1- 2*u + copula( u, u, ... ) ) - 1  # sec 3.3.2 Coles1999
  }, u=u, MoreArgs = list(...) )
  
  if( is.vector(barchiu)){
    plot(u, barchiu, xlab="u", ylab=expression(bar(chi)(u)))
  }  else {
    plot(0,0, xlab="u", ylab=expression(bar(chi)(u)),
         type = "n", xlim = range(u), ylim = range(barchiu))
    for(i in 1:nrow(barchiu)) lines(u, barchiu[i,])
  }
}


