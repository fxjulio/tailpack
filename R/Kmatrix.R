Kmatrix <- function( d=1, k , l ){
	npar <- k + l
	D <- matrix(0, npar-d, npar)
	for( i in 1:(npar-d) ){
		D[i,i+0:d] <- if(d==1) c(-1,1) else c(1,-2,1)
	}

	return( t(D)%*%D )
}
