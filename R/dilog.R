#'@title Dilogarithm function
#'
#'@description Computation of the dilogarithm function by nemerical integration.
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
