logistic <- function(x){
 x[ x > 16 ] <- 16
 v <- exp(x)/(1+exp(x))
 v[v>0.999] <- 0.999
 v[v<0.001] <- 0.001
 v
}

logit <- function(x) log(x/(1-x))