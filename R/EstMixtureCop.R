#'@title Estimation of bivariate mixture bivariate copula model
#'
#'@description Estimation of parameters from a mixture of  bivariate copula models
#'
#'@param y  (nx2) data matrix (observations or residuals) that will be transformed to pseudo-observations
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param reg   number of regimes
#'@param max_iter   maximum number of iterations of the EM algorithm
#'@param eps precision (stopping criteria); suggestion 0.0001.
#'
#'@author  Mamadou Yamar Thioub and Bruno Remillard, April 12, 2018
#'@return \item{theta}{(1 x reg) estimated parameter of the copula according to CRAN copula package (except for Frank copula, where theta = log(theta_R_Package)) for each component (except for degrees of freedom)}
#'@return \item{dof}{estimated degree of freedom, only for the Student copula}
#'@return \item{Q}{(1 x reg) estimated weights vector}
#'@return \item{eta}{(n x reg) conditional probabilities of being in regime k at time t given observations up to time t}
#'@return \item{tau}{estimated Kendall tau for each regime}
#'@return \item{U}{(n x 2) matrix of Rosenblatt transforms}
#'@return \item{cvm}{Cramer-von-Mises statistic for goodness-of-fit}
#'
#'@references <doi::10.1002/cjs.11534>
#'
#'
#'@export
EstMixtureCop<-function(y,reg,family,max_iter,eps){
  ninit=100
  n = dim(y)[1]; d = dim(y)[2]
  y = floor(  apply(y,2,rank) )/ (n + 1)
  r = reg*d
  n0 = floor((n/reg))
  ind0   = 1:n0
  alpha0 = rep(0,reg)
  tau    = rep(0,reg)

  for (j in 1:reg) {
    ind = (j-1)*n0 + ind0
    x   = y[ind, ]
    tau[j] = stats::cor(x,method = "kendall")[1,2]
    if (tau[j] <= 0.1){
      tau[j] = 0.1
    } else if (tau[j] >= 0.9) {
      tau[j] = 0.9
    }
  }
  alpha0 = ParamTau(family,tau)
  Q0 = rep(1,reg)/reg

  if (family == 't'){
    alpha0 = c(alpha0 , log(5))
  }

  # warm-up
  for (k in 1:ninit){
    emstep= EMstep(y,theta=alpha0,Q=Q0,family)
    alpha_new = emstep$theta_new; Qnew = emstep$Qnew; eta = emstep$eta
    Q0 = Qnew
    alpha0 = alpha_new
    #dof0 = dof_new
  }

  # iterations
  for (k in 1:max_iter){
    emstep= EMstep(y,alpha0,Q0,family)
    alpha_new = emstep$theta_new; Qnew = emstep$Qnew; eta = emstep$eta
    sum1 = sum(abs(alpha0))
    sum2 = sum(abs(alpha_new-alpha0))
    if ( sum2 < (sum1 * r * eps) ){
      break
    }
    Q0 = Qnew
    alpha0 = alpha_new
    #dof0 = dof_new;
  }

  # output
  alpha = alpha_new  # estimated unconstrained parameters
  Q = Qnew
  tau = KendallTau(family,alpha)
  dof = NaN
  theta = ParamCop(family,alpha)
  cdf = array(0,c(n,2,reg))
  switch(family,
         "gaussian" = {
           for( j in 1:reg){
             cdf[,,j] = RosenblattGaussian(y,theta[j])
           }
         },

         "t" = {
           dof = theta[(1+reg):length(theta)]
           for( j in 1:reg){
             cdf[,,j] = RosenblattStudent(y,theta[j],dof)
           }
         },

         "clayton" = {
           for( j in 1:reg){
             cdf[,,j] = RosenblattClayton(y,theta[j])
           }
         },

         "frank" = {
           for( j in 1:reg){
             cdf[,,j] = RosenblattFrank(y,exp(-theta[j]))
           }
         },

         "gumbel" = {
           for( j in 1:reg){
             cdf[,,j] = RosenblattGumbel(y,1/theta[j])
           }
         }
  )

  U=matrix(0,n,2)
  U[,1]= y[,1]
  for(j in 1:reg){
    U[,2] = U[,2]+ Q[j] * cdf[,2,j]
  }
  #V = floor(tiedrank(U)) / (n + 1)
  cvm = SnB(U)

  out = list(theta = theta, Q = Q, eta = eta, tau = tau, dof=dof, U = U, cvm = cvm)
  return(out)
}



EMstep <- function(y,theta,Q,family){
  n = dim(y)[1] #length of series
  r = length(Q) #number of regimes
  eta = matrix(0,n,r)
  c   = matrix(0,n,r)

  switch(family,
         "gaussian" = {
           for (j in 1:r){
             c[,j] = copulaFamiliesPDF('gaussian',y,theta[j])
           }
         },

         "t" = {
           for (j in 1:r){
             c[,j] = copulaFamiliesPDF('t',y,theta[j],theta[r+1])
           }
         },

         "clayton" = {
           for (j in 1:r){
             c[,j] = copulaFamiliesPDF('clayton',y,theta[j])
           }
         },

         "frank" = {
           for (j in 1:r){
             c[,j] = copulaFamiliesPDF('frank',y,theta[j])
           }
         },

         "gumbel" = {
           for (j in 1:r){
             c[,j] = copulaFamiliesPDF('gumbel',y,theta[j])
           }
         }
  )

  # eta
  v = c %*% diag(Q)
  v1= matrix(rep(rowSums(v),r), ncol=r)
  eta = v / v1
  Qnew = colMeans(eta)

  # finding theta
  switch(family,
         "gaussian" = {

           fun = function(thetaa){
             if (r < 2) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('gaussian',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('gaussian',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(eta[,i]* log( copulaFamiliesPDF('gaussian',y,thetaa[i]) )))
               }
               return(log_likelihood)
             }
           }

           if (r>= 2){
             theta_new =  stats::optim(par= theta,fun, method = "Nelder-Mead")$par
           } else if (r == 1){
             theta_new =  stats::optim(par= theta,fun, method = "L-BFGS-B")$par
           }
         },

         "t" = {

           fun = function(thetaa){
             if (r < 2) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('t',y,thetaa[1],thetaa[r+1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('t',y,thetaa[1],thetaa[r+1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(eta[,i]* log( copulaFamiliesPDF('t',y,thetaa[i],thetaa[r+1]) )))
               }
               return(log_likelihood)
             }
           }

           if (r>= 2){
             theta_new =  stats::optim(par= theta,fun, method = "Nelder-Mead")$par
           } else if (r == 1){
             theta_new =  stats::optim(par= theta,fun, method = "L-BFGS-B")$par
           }
         },

         "clayton" = {

           fun = function(thetaa){
             if (r < 2) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('clayton',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('clayton',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(eta[,i]* log( copulaFamiliesPDF('clayton',y,thetaa[i]) )))
               }
               return(log_likelihood)
             }
           }
           if (r>= 2){
             theta_new =  stats::optim(par= theta,fun, method = "Nelder-Mead")$par
           } else if (r == 1){
             theta_new =  stats::optim(par= theta,fun, method = "L-BFGS-B")$par
           }
         },

         "frank" = {

           fun = function(thetaa){
             if (r < 2) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('frank',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('frank',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(eta[,i]* log( copulaFamiliesPDF('frank',y,thetaa[i]) )))
               }
               return(log_likelihood)
             }
           }
           if (r>= 2){
             theta_new =  stats::optim(par= theta,fun, method = "Nelder-Mead")$par
           } else if (r == 1){
             theta_new =  stats::optim(par= theta,fun, method = "L-BFGS-B")$par
           }
         },

         "gumbel" = {

           fun = function(thetaa){
             if (r < 2) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('gumbel',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(eta[,1]* log( copulaFamiliesPDF('gumbel',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(eta[,i]* log( copulaFamiliesPDF('gumbel',y,thetaa[i]) )))
               }
               return(log_likelihood)
             }
           }
           if (r>= 2){
             theta_new =  stats::optim(par= theta,fun, method = "Nelder-Mead")$par
           } else if (r == 1){
             theta_new =  stats::optim(par= theta,fun, method = "L-BFGS-B")$par
           }
         }
  )
  out = list(theta_new=theta_new, Qnew=Qnew, eta=eta)
  return(out)
}
