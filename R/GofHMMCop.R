#'@title Goodness-of-fit of Markov regime switching bivariate copula model
#'
#'@description Goodness-of-fit test of a Markov regime switching bivariate copula model
#'
#'
#'@param       R         (n x 2) data matrix that will be transformed to pseudo-observations
#'@param       reg        number of regimes
#'@param       family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param       max_iter  maxmimum number of iterations of the EM algorithm
#'@param       eps       precision (stopping criteria); suggestion 0.0001
#'@param       n_sample  number of bootstrap; suggestion 1000
#'@param       n_cores    number of cores to use in the parallel computing
#'
#'@return \item{pvalue}{pvalue (significant when the result is greater than 5)}
#'@return \item{theta}{(1 x reg) estimated parameter of the copula according to CRAN copula package (except for Frank copula, where theta = log(theta_R_Package)) for each regime (except for degrees of freedom)}
#'@return \item{dof}{estimated degree of freedom, only for the Student copula}
#'@return \item{Q}{(reg x reg) estimated transition matrix}
#'@return \item{eta}{(n x reg) conditional probabilities of being in regime k at time t given observations up to time t}
#'@return \item{tau}{estimated Kendall tau for each regime}
#'@return \item{U}{(n x 2) matrix of Rosenblatt transforms}
#'@return \item{cvm}{Cramer-von-Mises statistic for goodness-of-fit}
#'@return \item{W}{regime probabilities for the conditional distribution given the past Kendall's tau}
#'
#'
#'@references <doi::10.1002/cjs.11534>
#'
#@examples Q <- matrix(c(0.8, 0.2, 0.3, 0.7),2,2) ; kendallTau <- c(0.3 ,0.7) ;
#data <- SimHMMCop(Q, 'clayton', kendallTau, 10)$SimData;
#gof <- GofHMMCop(data,2,'clayton',10000,0.0001,1)
#'
#'
#'@export
#'
GofHMMCop <-function(R, reg, family, max_iter ,eps ,n_sample,n_cores){
  registerDoParallel(n_cores)

  esthmcop =  EstHMMCop(R, reg, family, max_iter, eps)
     theta = esthmcop$theta
         Q = esthmcop$Q
       eta = esthmcop$eta
        nu = esthmcop$nu
       tau = esthmcop$tau
        df = esthmcop$dof
         U = esthmcop$U
   cvm_est = esthmcop$cvm
         W = esthmcop$W
  eta0 = sample(1:reg, n_sample, replace = T)
  n = dim(R)[1];


  #result <- foreach(i=1:n_sample) %do% bootstrapfun(Q,family,tau,n,df,max_iter,eps,1)
  result <- foreach(i=1:n_sample, .packages='HMMcopula') %dopar% bootstrapfun(Q,family,tau,n,df,max_iter,eps,1)

  cvm_sim1 = rep(0,n_sample)
  for (i in 1:n_sample){
    cvm_sim1[i] = result[[i]]$cvm_sim
  }

  pvalue = 100*mean( cvm_sim1 > cvm_est)

  out = list(  pvalue =  pvalue, theta = theta, Q = Q, eta = eta, nu = nu, tau = tau, df = df, U = U, cvm_est = cvm_est, W = W)
  return(out)
}
