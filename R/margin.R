margin <- function( data, t=seq(1,nrow(data)), logret=FALSE, tail=-1, mar="unif", filtered=FALSE){
	if( logret ){
    lret <- apply(data, 2, function(x) diff(log(x)) )
    out.t <- t[-1]
  } else {
    lret <- data
    out.t <- t
  }

	flag <- apply(lret, 1, function(x) !any(is.na(x)))
	lret <- lret[flag,]
	out.t <- out.t[flag]

	flag <- apply( lret, 1, function(x) all( x != 0 ) )
	lret <- lret[flag,]
	out.t <- out.t[flag]

	if( filtered ){
		lret <- apply( lret, 2, function(x){
			fit <- garch(x, order = c(1, 1), control=garch.control(trace = FALSE))
			#summary(fit)
			residuals(fit)
		})
		lret <- lret[-1,]
		out.t <- out.t[-1]
	}

	#flag <- apply( lret, 1, function(x) all( tail*x > 0 ))
	#tret <- tail*lret[flag,]
	#out.t <- out.t[flag]

	tret <- tail*lret

	n <- nrow(tret)+1
	uret <- apply( tret, 2, rank)/n 

	attr(uret, 't') <- out.t
	attr(tret, 't') <- out.t

	if( mar=="" ){
		return( tret )
	}
	if( mar=="unif" )
		return( uret )
	if( mar=="frechet" )
		return( -1/log(uret) )
	if( mar=="pareto" )
		return( 1/(1-uret) )
	if( mar=="exponential" )
		return( -log( 1-uret ) )

}
