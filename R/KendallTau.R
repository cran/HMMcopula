#'@title Kendall's tau of a copula
#'
#'@description This function computes the Kendall's tau of a copula family with a unconstrainted parameter alpha.
#'
#'@param family "gaussian" , "t" , "clayton" , "frank" , "gumbel"
#'@param alpha  unconstrainted parameters of the copula family
#'
#'@return \item{tau}{estimated Kendall's tau}
#'
#'@export
KendallTau<- function(family,alpha){

  switch(family,
         "gaussian" = { Rho = ( 2*exp(alpha) / (exp(alpha)+1) ) - 1
         Rho[is.na(Rho)]=1
         tau = (2 / pi) * asin(Rho)
         },

         "t" = {
           Rho = ( 2*exp(alpha[1:(length(alpha)-1)]) / (exp(alpha[1:(length(alpha)-1)])+1) ) - 1
           Rho[is.na(Rho)]=1
           tau = (2 / pi) * asin(Rho)
         },

         "clayton" = {
           theta = exp(alpha)
           tau = theta / (theta+2)
         },

         "frank" = {
           theta = exp(-alpha)

           tau = (log(theta)^2 + 4*log(theta) + 4*sapply(theta,dilog)) / (log(theta)^ 2)
         },

         "gumbel" = {
           theta = 1 / (exp(alpha)+1)
           tau = 1 - theta
         }
  )

  return(tau)
}


# dilog = function(x){
#   if(x==1) { return(NA)}
#   out = stats::integrate(function(y) log(y)/(1-y), lower = 1, upper = x)$value
#   return(out)
# }







