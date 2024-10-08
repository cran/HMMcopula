#'@title Simulation of bivariate Markov regime switching copula model
#'
#'@description Simulaion of values from a bivariate Markov regime switching copula model
#'
#'@param Q  Transition probality matrix (d x d);
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param KendallTau   Kendall's rank correlation
#'@param n    number of simulated vectors
#'@param DoF degree of freedom only for the Student copula
#'
#'@return \item{SimData}{Simulated Data}
#'@return \item{MC}{Markov chain regimes}
#'@return \item{alpha}{parameters alpha}
#'
#'@examples Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) ; kendallTau <- c(0.3 ,0.7) ;
#'simulations <- SimHMMCop(Q, 'gumbel', kendallTau, 300)
#'
#'
#'@export
SimHMMCop<-function(Q, family, KendallTau, n, DoF){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }

  #if(dim(Q)[2] <= 1){
    # warning("Transition matrix not assigned, will simulate one regime copula.")
  #}else
  if(nargs()<=3 ){
    stop("Requires at least four input arguments.")
  } else if(nargs()<=4 & family == 't'){
    stop("Requires at five input arguments.")
  }




  if(reg >=2){
    MC = SimMarkovChain(Q,n)
  }else  MC = rep(1,n+1)
  alpha = rep(0,reg)
  Sim   = matrix(0,n,2*reg)
  SimData = matrix(0,n,2)

  for(k in 1:reg){

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

           for( k in 1:reg){
             u = copula::rCopula(n,copula::normalCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }

         },

         "t" = {
           for(k in 1:reg){
             #if (regd == 1){
               u = copula::rCopula(n, copula::tCopula(alpha[k], dim = 2,df = DoF))
             #} else {
              # u = copula::rCopula(n, copula::tCopula(alpha[k], dim = 2,df = DoF[k]))
             #}
             Sim[,(2*k-1):(2*k)] = u
           }
         },

         "clayton" = {
           for( k in 1:reg){
             u = copula::rCopula(n,copula::claytonCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }
         },

         "frank" = {
           for( k in 1:reg){
             u = copula::rCopula(n,copula::frankCopula(alpha[k], dim = 2))
             Sim[,(2*k-1):(2*k)] = u
           }
         },

         "gumbel" = {
           for( k in 1:reg){
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

