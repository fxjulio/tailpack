calcKnots <- function( x, k=20, l=3){
	h <- (max(x)-min(x))/k
	
	seq( min(x)-l*h, max(x)+l*h, by=h )

}
