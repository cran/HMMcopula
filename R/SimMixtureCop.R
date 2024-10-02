#'@title Simulation of bivariate mixture copula model
#'
#'@description Simulation of observations from a bivariate mixture copula model
#'
#'@param Q Weights vector (1 x component);
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param  KendallTau   Kendall's rank correlation
#'@param  n    number of simulated vectors
#'@param DoF vector of degree of freedom only for the Student copula
#'
#'@return \item{SimData}{Simulated Data}
#'@return \item{MC}{Markov chain regimes}
#'@return \item{alpha}{parameters alpha}
#'
#'@examples Q <- matrix(c(0.8, 0.2),1,2) ; kendallTau <- c(0.3 ,0.7) ;
#'simulations <- SimMixtureCop(Q, 'gaussian', kendallTau, 300)
#'
#'
#'
#'@export
#'
SimMixtureCop<-function(Q, family, KendallTau, n, DoF){

if(is.null(dim(Q))){
  QQ0 = matrix(Q)
  reg = dim(QQ0)[1]
} else {
  reg = dim(Q)[2]
}
# if(is.null(dim(DoF))){
#   DoF0 = matrix(DoF)
#   regd = dim(DoF0)[1]
# } else {
#   regd = dim(DoF)[2]
# }

if(nargs()<=3){ stop('Requires at least four input arguments.')
}else if(round(max(sum(Q))) != 1 | round(min(sum(Q))) != 1){
    stop('Wrong transition matrix.')
  }

# number of regimes in the simulated data.

if(reg >= 2){
  MC =  sample(1:reg, n,TRUE, t(Q)) #simulated regimes

} else {
  MC = matrix(0,n+1)
  MC[1:(n+1),1] = 1 } # MC[1:(n+1),1] = 1  # in the case of one regime copula

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
