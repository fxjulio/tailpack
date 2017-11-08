#' Hill estimator
#' @param x Vector of values
#' @param k Number of order statistics to consider
#' @export
hill <- function(x, k){
  sx <- sort(x, decreasing=TRUE )
  mean( log(sx[1:k]/sx[k+1])  )
}
