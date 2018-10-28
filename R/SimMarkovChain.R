#'Markov chain simulation
#'
#'@description This function generates a Markov chain X(1), ..., X(n) with transition matrix Q,
#' starting from a state eta0 or the uniform distribution on {1,..., r}
#'
#'@param  Q Transition probality matrix (d x d)
#'@param  n number of simulated vectors
#'@param  eta0 variable eta
#'
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

