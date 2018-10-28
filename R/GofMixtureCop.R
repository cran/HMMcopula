#'@title Goodness-of-fit of mixture bivariate copula model
#'
#'@description This function performs goodness-of-fit test of a mixture bivariate copula model
#'
#' @param       R         (nx2) data matrix (observations or residuals) that will be transformed to pseudo-observations
#' @param       reg        number of regimes
#' @param       family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#' @param       max_iter    maxmimum number of iterations of the EM algorithm
#' @param       eps       precision (stopping criteria); suggestion 0.0001
#' @param       n_sample    number of bootstrap; suggestion 1000
#' @param       n_cores    number of cores to use in the parallel computing
#'
#'@return \item{pvalue}{pvalue (significant when the result is greater than 5)}
#'@return \item{theta}{(1 x reg) estimated parameter of the copula according to CRAN copula package (except for Frank copula, where theta = log(theta_R_Package)) for each component (except for degrees of freedom)}
#'@return \item{dof}{estimated degree of freedom, only for the Student copula}
#'@return \item{Q}{(1 x reg) estimated weights vector}
#'@return \item{eta}{(n x reg) conditional probabilities of being in regime k at time t given observations up to time t}
#'@return \item{tau}{estimated Kendall tau for each regime}
#'@return \item{U}{(n x 2) matrix of Rosenblatt transforms}
#'@return \item{cvm}{Cramer-von-Mises statistic for goodness-of-fit}
# @author  By Bruno Remillard, Nov 28, 2010
#'
#'
#'
#'@export
GofMixtureCop<-function(R,reg,family,max_iter,eps,n_sample,n_cores){
  registerDoParallel(n_cores)

  esthmmcop =  EstMixtureCop(R,reg,family,max_iter,eps)

  theta=esthmmcop$theta; Q=esthmmcop$Q; eta=esthmmcop$eta;
  tau=esthmmcop$tau; df=esthmmcop$dof; U=esthmmcop$U; cvm_est=esthmmcop$cvm

  #cvm_sim = rep(0,n_sample)
  eta0 = sample(1:reg,n_sample,replace=T)

  n = dim(R)[1];

  # parametric bootstrap


  result <- foreach(i=1:n_sample, .packages='HMMcopula') %dopar% bootstrapfun(Q,family,tau,n,df,max_iter,eps,0)


  cvm_sim1 = rep(0,n_sample)
  for (i in 1:n_sample){
    cvm_sim1[i] = result[[i]]$cvm_sim
  }

  #for(i in 1:n_sample){
  #  R1 = SimMixtureCop(Q, family, tau, n, df)$SimData
  #  esthmmco = EstMixtureCop(R1,reg,family,max_iter,eps)
  #  theta1=esthmmco$theta ; Q1=esthmmco$Q ; eta1=esthmmco$eta ; tau1=esthmmco$tau ; dof1=esthmmco$df ; Usim=esthmmco$U ; cvm_sim[i]=esthmmco$cvm
  #}

  #output
  pvalue = 100*mean( cvm_sim1 > cvm_est)

  out = list( pvalue=  pvalue, theta=theta, Q=Q, eta=eta, tau=tau, df=df, U=U, cvm_est=cvm_est)
  return(out)
}
