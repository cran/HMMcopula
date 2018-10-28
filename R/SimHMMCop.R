#'@title Simulation of bivariate Markov regime switching copula model
#'
#' @description This function simulates observation from a bivariate Markov regime switching copula model
#'
#'@param Q  Transition probality matrix (d x d);
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param KendallTau   Kendall's rank correlation
#'@param n    number of simulated vectors
#'@param DoF degree of freedom only for the Student copula
#@param DoF vector of degree of freedom (1 x d) or one degree of freedom (for all regime),
#'
#'@return \item{SimData}{Simulated Data}
#'@return \item{MC}{Markov chain regimes}
#'@return \item{alpha}{parameters alpha}
#'
#'@examples Q <- matrix(c(0.8, 0.2, 0.3, 0.7),2,2) ; kendallTau <- c(0.3 ,0.7) ;
#'simulations <- SimHMMCop(Q, 'gumbel', kendallTau, 300)
#'
#'
#'@export
SimHMMCop<-function(Q, family, KendallTau, n, DoF){
  if(dim(Q)[2] <= 1){
    warning("Transition matrix not assigned, will simulate one regime copula.")
  }else if(nargs()<=3 ){
    stop("Requires at least four input arguments.")
  } else if(nargs()<=4 & family == 't'){
    stop("Requires at five input arguments.")
  }

  # if(is.null(dim(DoF))){
  #   DoF0 = matrix(DoF)
  #   regd = dim(DoF0)[1]
  # } else {
  #   regd = dim(DoF)[2]
  # }

  Nregimes = dim(Q)[2]
  if(Nregimes >=2){
    MC = SimMarkovChain(Q,n)
  }else  MC = rep(1,n+1)
  alpha = rep(0,Nregimes)
  Sim   = matrix(0,n,2*Nregimes)
  SimData = matrix(0,n,2)

  for(k in 1:Nregimes){

    switch(family,
           "gaussian" = {
             alpha[k] =  copula::iTau(copula::normalCopula(),  tau = KendallTau[k])
           },

           "t" = {
             alpha[k] =  copula::iTau(copula::normalCopula(),  tau = KendallTau[k]) #copula::iTau(copula::tCopula(df = DoF ),  tau = KendallTau[k])
           },

           "clayton" = {
             alpha[k] =  copula::iTau(copula::claytonCopula(),  tau = KendallTau[k])
           },

           "frank" = {
             alpha[k] =  copula::iTau(copula::frankCopula(),  tau = KendallTau[k])
           },

           "gumbel" = {
             alpha[k] =  copula::iTau(copula::gumbelCopula(),  tau = KendallTau[k])
           }
    )



  }

  switch(family,
         "gaussian" = {

           for( k in 1:Nregimes){
             u = copula::rCopula(n,copula::normalCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }

         },

         "t" = {
           for(k in 1:Nregimes){
             #if (regd == 1){
               u = copula::rCopula(n, copula::tCopula(alpha[k], dim = 2,df = DoF))
             #} else {
              # u = copula::rCopula(n, copula::tCopula(alpha[k], dim = 2,df = DoF[k]))
             #}
             Sim[,(2*k-1):(2*k)] = u
           }
         },

         "clayton" = {
           for( k in 1:Nregimes){
             u = copula::rCopula(n,copula::claytonCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }
         },

         "frank" = {
           for( k in 1:Nregimes){
             u = copula::rCopula(n,copula::frankCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }
         },

         "gumbel" = {
           for( k in 1:Nregimes){
             u = copula::rCopula(n,copula::gumbelCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }
         }
  )


  for(i in 1:n){
    k = MC[i]
    SimData[i,] = Sim[i,(2*k-1):(2*k)]
  }
  out=list(SimData=SimData,MC=MC,alpha=alpha)
  return(out)
}

