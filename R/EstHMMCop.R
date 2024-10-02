#'@title Estimation of bivariate Markov regime switching bivariate copula model
#'
#'@description Estimation of parameters from a bivariate Markov regime switching bivariate copula model
#'
#'@param y   (nx2) data matrix (observations or residuals) that will be transformed to pseudo-observations
#'@param family    'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'
#'@param reg    number of regimes
#'@param max_iter  maximum number of iterations of the EM algorithm
#'@param eps precision (stopping criteria); suggestion 0.0001.
#'
#'@author Mamadou Yamar Thioub and Bruno Remillard, April 12, 2018
#'@return \item{theta}{(1 x reg) estimated parameter of the copula according to CRAN copula package (except for Frank copula, where theta = log(theta_R_Package)) for each regime (except for degrees of freedom)}
#'@return \item{dof}{estimated degree of freedom, only for the Student copula}
#'@return \item{Q}{(reg x reg) estimated transition matrix}
#'@return \item{eta}{(n x reg) conditional probabilities of being in regime k at time t given observations up to time t}
#'@return \item{tau}{estimated Kendall tau for each regime}
#'@return \item{U}{(n x 2) matrix of Rosenblatt transforms}
#'@return \item{cvm}{Cramer-von-Mises statistic for goodness-of-fit}
#'@return \item{W}{regime probabilities for the conditional distribution given the past Kendall's tau}
#'
#'@references <doi::10.1002/cjs.11534>
#'
#'
#'@examples Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) ; kendallTau <- c(0.3 ,0.7) ;
#'data <- SimHMMCop(Q, 'clayton', kendallTau, 10)$SimData;
#'estimations <- EstHMMCop(data,2,'clayton',10000,0.0001)
#'
#'
#'@export
EstHMMCop<-function(y,reg,family,max_iter,eps){
  ninit=100   #minimum number of iterations
  n = dim(y)[1]; d = dim(y)[2]
  y = floor( apply(y,2,rank) )/ (n + 1)
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
  Q0 = matrix(1,reg,reg)/reg

  if (family == 't'){
    alpha0 = c(alpha0 , log(5))
  }

  # warm-up
  for (k in 1:ninit){
    emstep= EMStep(y=y,theta=alpha0,Q=Q0,family=family)
    nu=emstep$nu ;alpha_new = emstep$theta_new; Qnew = emstep$Qnew; eta = emstep$eta
    eta_bar =emstep$eta_bar; lambda= emstep$lambda; Lambda=emstep$Lambda
    Q0 = Qnew
    alpha0 = alpha_new
    #dof0 = dof_new
  }

  # iterations
  for (k in 1:max_iter){
    emstep= EMStep(y=y,theta=alpha0,Q=Q0,family=family)
    nu = emstep$nu ;alpha_new = emstep$theta_new; Qnew = emstep$Qnew; eta = emstep$eta
    eta_bar = emstep$eta_bar; lambda= emstep$lambda; Lambda = emstep$Lambda

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
           dof = theta[(1+reg)]
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
  eta00 = rep(1,reg)/reg

  w00 = rbind(eta00,eta) %*% Q

  W   = w00[1:(dim(w00)[1]-1),]
  if (reg == 1){
    W = data.frame(as.list(W))
  }

  U=matrix(0,n,2)
  U[,1]= y[,1]
  for(j in 1:reg){
    U[,2] = U[,2]+ W[,j] * cdf[,2,j]
  }
  #V = floor(tiedrank(U)) / (n + 1)
  cvm = SnB(U)
  out = list(theta=theta, Q=Q, eta=eta, nu=nu,tau=tau, dof=dof, U=U, cvm=cvm,W=W)
  return(out)
}



EMStep <- function(y,theta,Q,family){
  n = dim(y)[1] #length of series
  r = dim(Q)[2] #number of regimes
  eta_bar = matrix(0,n,r)
  eta = matrix(0,n,r)
  lambda   = matrix(0,n,r)
  c   = matrix(0,n,r)
  Lambda   = array(0,c(r,r,n))
  M        = matrix(0,r,r)
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

  # eta_bar
  eta_bar[n,]=1/r

  for(k in 1:(n-1)){
    i = n-k        # backward
    j = i+1        # j is the index of period (t+1)
    v =  ( eta_bar[j,] * c[j,] ) %*% t(Q) #  numerator of eta_bar(i)
    eta_bar[i,] = v/sum(v);
  }


  # eta
  eta0 = rep(1,r)/r

  v = ( eta0 %*% Q) * c[1,]  #numerator of eta at  t = 1;
  eta[1,] = v/sum(v)

  for (i in 2:n){
    v = ( eta[i-1,] %*% Q) * c[i,]    # numerateur du eta lorsque t > 1;
    eta[i,] = v/sum(v)        #valeur du eta lorsque t > 1;

  }
  #lambda
  v      = eta * eta_bar    # numerator of lambda
  sv0    = rowSums(v)       # calcul des denomirateur des lambda

  for(j in 1:r){
    lambda[,j] = v[,j] / sv0        # calcul des lambda
  }

  #lambda
  gc =  eta_bar * c   #les deux derniers calcul (plus \E0 droite) du numerateur de Lambda

  M = Q * (as.matrix(eta0)%*% gc[1,] ) # numerator of Lambda
  MM = sum(M)   #denominator of Lambda
  Lambda[,,1] = M/MM

  for(i in 2:n){
    M = Q * ( as.matrix(eta[i-1,]) %*% gc[i,]  )    # numerateur de Lambda
    MM = sum(M)    #denominateur de Lambda
    Lambda[,,i] = M/MM    #Lambda pour t < n
  }
  nu = colMeans(lambda)
  Qnew = Q

  if (r >= 2){
    #for(j in 1:r){
    #Lambda2_1 = matrix(Lambda[j,1,], n, 1)
    #Lambda2_2 = matrix(Lambda[j,2,], n, 1)
    #Lambda2 = matrix(c(Lambda2_1,Lambda2_2),n,2)
    #sv = apply(Lambda2, 2, sum)
    #####sv = rowSums(Lambda2, dims=2)
    ####sv = rowSums(Lambda[j,,])
    #ssv = sum(sv)
    #Qnew[j,] = sv / ssv
    #}
    for(j in 1:r){
      sv = rowSums(Lambda[j,,], dims=1)
      ssv = sum(sv)
      Qnew[j,] = sv/ssv
    }
  }




  # finding theta
  switch(family,
         "gaussian" = {

           fun = function(thetaa){
             if (r < 2) {
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('gaussian',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('gaussian',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(lambda[,i]* log( copulaFamiliesPDF('gaussian',y,thetaa[i]) )))
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
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('t',y,thetaa[1],thetaa[r+1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('t',y,thetaa[1],thetaa[r+1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(lambda[,i]* log( copulaFamiliesPDF('t',y,thetaa[i],thetaa[r+1]) )))
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
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('clayton',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('clayton',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(lambda[,i]* log( copulaFamiliesPDF('clayton',y,thetaa[i]) )))
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
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('frank',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('frank',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(lambda[,i]* log( copulaFamiliesPDF('frank',y,thetaa[i]) )))
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
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('gumbel',y,thetaa[1]) ))
             } else if (r > 1) {
               log_likelihood =  -sum(lambda[,1]* log( copulaFamiliesPDF('gumbel',y,thetaa[1]) ))
               for(i in 2:r){
                 log_likelihood =log_likelihood + (-sum(lambda[,i]* log( copulaFamiliesPDF('gumbel',y,thetaa[i]) )))
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
  out = list(nu=nu, theta_new=theta_new, Qnew=Qnew, eta=eta, eta_bar=eta_bar, lambda=lambda, Lambda=Lambda)
  return(out)
}
