

etalogistic1 <- function(n, alpha, eta, rho=1){

	U <- matrix(runif(2*n), n, 2)

	Nrho <- rho^(-1/eta) + rho^(1/eta) - (rho^(-1/alpha)+rho^(1/alpha))^(alpha/eta)
	uast <- Nrho^(-1)*( rho^(-1/eta) - rho^(-1/alpha)*( rho^(-1/alpha) + rho^(1/alpha))^(alpha/eta-1) )
	
	flag <- U[,1] <= uast

	V<-numeric(n)
	V[flag] <- acos( (1-Nrho*rho^(1/eta)*U[flag,1])^(eta/(2*(eta-alpha)) ) )
	V[!flag] <- asin( (rho^(-1/eta)*(Nrho*U[!flag,1]+(rho^(-1/alpha)+rho^(1/alpha))^(alpha/eta)-rho^(-1/eta)))^(eta/(2*(eta-alpha))) )

	Z<-numeric(n)
	Z[flag] <- U[flag,2]*rho^(-1/eta)*cos(V[flag])^(-2*alpha/eta)
	Z[!flag] <- U[!flag,2]*rho^(1/eta)*sin(V[!flag])^(-2*alpha/eta)

	#plot(V,Z)

	cbind( S=rho^(-1)*Z^(-eta)*cos(V)^(-2*alpha), T=rho*Z^(-eta)*sin(V)^(-2*alpha) )

}

etalogistic2 <- function(n, alpha, eta, rho=1){

	Z <- replicate(n, ifelse( runif(1) < alpha/eta, rgamma(1,2,1), rexp(1,1) ) )
	V <- asin( sqrt(runif(n)) )


	Nrho <- rho^(-1/eta) + rho^(1/eta) - (rho^(-1/alpha)+rho^(1/alpha))^(alpha/eta)

	# x,y >0 block maximum
	cbind( X=rho^(-1)*(Nrho*Z)^(-eta)*cos(V)^(-2*alpha), Y=rho*(Nrho*Z)^(-eta)*sin(V)^(-2*alpha) )

}