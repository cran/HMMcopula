#'@title Markov chain simulation
#'
#'@description Simulation of  n consecutive  values of a Markov chain with transition matrix Q, 
#' starting from a state eta0 or the uniform distribution on the set 1,..., r.
#'
#'@param  Q Transition probability matrix (d x d)
#'@param  n number of simulated vectors
#'@param  eta0 variable eta
#'
#'@return \item{x}{Simulated Markov chain sequence}
#'@export
#'
SimMarkovChain<-function(Q,n,eta0){
  r = dim(Q)[1] ; p = dim(Q)[2]
  x = rep(0,n)  ; x0= matrix(0,n,r)
  if(nargs()<3){
    ind = sample(1:r,1)
  }else ind = eta0

  for( k in 1:r){
    x0[,k] = sample(1:r, size = n, replace = TRUE, prob = Q[k,])
  }

  for( i in 1:n){
    x[i] = x0[i,ind]
    ind = x[i]
  }
  return(x)
}

