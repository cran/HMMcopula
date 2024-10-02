#'@title Bootstrap function for a  bivariate copula models
#'@description Bootstrapping function needed for parallel computing
#'
#'@param Q Weights vector (1 x reg or component);
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param tau   Kendall's rank correlation
#'@param n    number of simulated vectors
#'@param df vector of degree of freedom (d x 1), only for the Student copula.
#'@param max_iter maximum number of iterations for estimation
#'@param eps       precision (e.g 0.00001);
#'@param HMM    1 (if HMM) , 0 (if mixture);
#'
#'@author  Mamadou Yamar Thioub and Bruno Remillard, April 12, 2018
#'@return \item{theta1}{Estimated copula parameters}
#'@return \item{Q1}{Estimated transition matrix}
#'@return \item{eta1}{Estimated probabilites for regimes}
#'@return \item{tau1}{Estimated Kendall's tau}
#'@return \item{dof1}{Estimated degrees of freedom for the Student copula}
#'@return \item{Usim}{Estimated pseudo-observations}
#'@return \item{cvm_sim}{Estimated Cramer-von Mises statistic}
#'@export
#'@keywords internal

bootstrapfun <- function(Q,family,tau,n,df,max_iter,eps,HMM){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }

  if (HMM == 0){R1 = SimMixtureCop(Q, family, tau, n, df)$SimData
    esthmmco = EstMixtureCop(R1,reg,family,max_iter,eps)
    out = list(theta1=esthmmco$theta , Q1=esthmmco$Q , eta1=esthmmco$eta , tau1=esthmmco$tau , dof1=esthmmco$df , Usim=esthmmco$U , cvm_sim=esthmmco$cvm)
    return(out)

  } else if (HMM == 1) {R1 = SimHMMCop(Q, family, tau, n, df)$SimData
    esthmmco = EstHMMCop(R1, reg, family, max_iter, eps)
    out = list(theta1=esthmmco$theta , Q1=esthmmco$Q , eta1=esthmmco$eta , nu1 = esthmmco$nu , tau1=esthmmco$tau , dof1=esthmmco$df , Usim=esthmmco$U , cvm_sim=esthmmco$cvm)
    return(out)

  }

}
