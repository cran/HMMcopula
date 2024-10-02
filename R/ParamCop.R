#'@title Theta estimation
#'
#'@description Parameters of a copula according to CRAN copula package (except for Frank copula, where theta = log(theta_R_Package)), corresponding to the unconstrainted parameters alpha.
#'
#'@param family "gaussian" , "t" , "clayton" , "frank" , "gumbel"
#'@param alpha  unconstrainted parameters of the copula family
#'
#'@return \item{theta}{matlab parameters}
#'
#'@export
ParamCop<- function(family,alpha){

  switch(family,
         "gaussian" = {
           theta =  2 /(1+exp(-alpha) )-1
         },

         "t" = {
           reg = length(alpha) - 1
           theta = 2 /(1+exp(-alpha[1:reg])) -1
           dof   = exp(alpha[reg+1])
           theta = c(theta,dof)
         },

         "clayton" = {
           theta = exp(alpha)
         },

         "frank" = {
           theta = alpha
         },

         "gumbel" = {
           theta = exp(alpha)+1
         }
  )

  return(theta)
}





