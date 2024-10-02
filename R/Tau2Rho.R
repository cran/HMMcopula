#'@title Spearman's rho
#'
#'@description Value  of Spearman's rho corresponding to a constrainted (matlab) parameter theta for a copula family.
#'
#'@param family  'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param   theta  parameter of the copula according to CRAN copula package (except for Frank copula, where theta = log(theta_R_Package))
#'
#'@return  \item{rho}{estimated Spearman's rho}
#'
#'
#'@export
#'
Tau2Rho <-function(family,theta){
  N=100000
  # if (require("copula", quietly=TRUE)) {
  #   copula::iTau
  #   copula::normalCopula
  #   copula::claytonCopula
  #   copula::frankCopula
  #   copula::gumbelCopula
  #   copula::tCopula
  #   copula::rCopula
  # } else {
  #     stop("Please install package 'pkgB' to do this.")
  # }



  switch(family,
         "gaussian" = {
           u =   copula::rCopula(N,copula::normalCopula(theta[1], dim = 2))

         },

         "t" = {
           u =   copula::rCopula(N,copula::tCopula(theta[1], df=theta[2], dim = 2))
         },

         "clayton" = {
           u =   copula::rCopula(N,copula::claytonCopula(theta[1], dim = 2))
         },

         "frank" = {
           u =   copula::rCopula(N,copula::frankCopula(theta[1], dim = 2))
         },

         "gumbel" = {
           u =   copula::rCopula(N,copula::gumbelCopula(theta[1], dim = 2))
         }
  )

  rho  <- stats::cor(u, method = "spearman")[1,2]
  return(rho)
}

