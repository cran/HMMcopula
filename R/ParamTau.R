#'@title Alpha estimation
#'
#'@description Unconstrainted parameter for a given Kendall's tau.
#'
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param tau        Kendall's tau of the copula family
#'
#'@return \item{alpha}{estimated unconstrainted parameter}
#'
#'@export
ParamTau<- function(family,tau){
  reg = length(tau)
  alpha = rep(0,reg)
  for (j in 1:reg){
    switch(family,
           "gaussian" = {
             alpha[j] =  copula::iTau(copula::normalCopula(),  tau = tau[j])
           },

           "t" = {
             alpha[j] =  copula::iTau(copula::normalCopula(),  tau = tau[j]) #iTau(tCopula(df = df ),  tau = tau[j])
           },

           "clayton" = {
             alpha[j]=  copula::iTau(copula::claytonCopula(),  tau = tau[j])
           },

           "frank" = {
             alpha[j] =  copula::iTau(copula::frankCopula(),  tau = tau[j])
           },

           "gumbel" = {
             alpha[j] =  copula::iTau(copula::gumbelCopula(),  tau = tau[j])
           }
    )

  }



  switch(family,
         "gaussian" = {
           alpha = log( (1+alpha) / (1-alpha) )
         },

         "t" = {
           alpha = log( (1+alpha) / (1-alpha) )
         },

         "clayton" = {
           alpha = log(alpha)
         },

         "frank" = {
           alpha = alpha
         },

         "gumbel" = {
           alpha = log(alpha-1)
         }
  )

  return(alpha)
}
