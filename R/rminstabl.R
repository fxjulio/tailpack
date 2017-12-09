rminstabl <- function(n, x){
  covmat <- exp( -as.matrix(dist(x, method = "manhattan"))  )
  
  W <- replicate(n, {
    W <- rep(1e20, length(x))
    poisson <- 0
    nOK <- length(x)
    
    while (nOK){
      poisson <- poisson + rexp(1)
      thresh <- poisson/uBound
      
      Y <- rmvnorm(1, rep(0, length(x)), covmat)
      W <- pmin( W,  poisson/pmax(1e-10,Y) )
      
      nOK <- length(x) - sum(thresh >= W)
    }
    W
  })
  
  1/sqrt(2*pi)*t(W)
}