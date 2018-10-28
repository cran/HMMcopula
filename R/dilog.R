#'@title Dilogarithm function
#'
#'@description This function computes the dilogarithm of a number.
#'
#'@param x a real number
#'
#'@return \item{out}{dilogarithm}
#'
#'@export
dilog <- function(x){
  if(x==1) { return(NA)}
  out = stats::integrate(function(y) log(y)/(1-y), lower = 1, upper = x)$value
  return(out)
}
