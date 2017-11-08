#' Logistic samples
#' @export
rbvlog_shi <- function(n, alpha){
  sim <- matrix(0,n,2)
  u <- runif(n)
  z <- rexp(n)
  flag <- runif(n) < alpha
  z[flag] <- rexp(sum(flag))+rexp(sum(flag))

  sim[,1] <- 1/(z * u^alpha)
  sim[,2] <- 1/(z * (1-u)^alpha)

  sim

}
