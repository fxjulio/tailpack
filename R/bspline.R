bspline <- function( s, xi, i, l ){
	h <- xi[2]-xi[1]
	if( l==0 ) as.numeric(  xi[i]<= s & s < xi[i+1]  )
	else {
		(s-xi[i])/(l*h)*bspline( s, xi, i, l-1) + 
		(xi[i+l+1]-s)/(l*h)*bspline( s, xi, i+1, l-1 )
	}
}


Bmat <- function( x, xi, k, l ){
	X <- matrix(0, length(x), k+l)
	for( i in 1:ncol(X) ) X[,i] <- bspline( x, xi, i, l)
	X
}

bsIntercept <- function( b, k, l){
	if( l == 2 ){
		w <- rep(1, k+l)
		w[1] <- w[k+l] <- 1/6
		w[2] <- w[k+l-1] <- 5/6
	}

	if( l == 3 ){
		w <- rep(1, k+l)
		w[1] <- w[k+l] <- 1/24
		w[2] <- w[k+l-1] <- 12/24
		w[3] <- w[k+l-2] <- 23/24
	}
  
  return( sum(b*w)/k )	

}
