#'@title CopulaFamiliesCDF
#'
#'@description COPULACDF Cumulative probability function for a copula with linear correlation parameters RHO
#'
#' @param family copula familly= "gaussian" , "t" , "clayton" , "frank" , "gumbel"
#' @param u is an N-by-P matrix of values in [0,1], representing N points in the P-dimensional unit hypercube
#' @param ...  additionnal parameter like RHO  a P-by-P correlation matrix.
#'
#' @return
#'Y = COPULACDF('Gaussian',U,RHO) returns the cumulative probability of the
#'   Gaussian copula with linear correlation parameters RHO, evaluated at the
#'   points in U. U is an N-by-P matrix of values in [0,1], representing N
#'   points in the P-dimensional unit hypercube.  RHO is a P-by-P correlation
#'   matrix.  If U is an N-by-2 matrix, RHO may also be a scalar correlation
#'   coefficient.
#'
#'   Y = COPULACDF('t',U,RHO,NU) returns the cumulative probability of the t
#'   copula with linear correlation parameters RHO and degrees of freedom
#'   parameter NU, evaluated at the points in U.  U is an N-by-P matrix of
#'   values in [0,1]. RHO is a P-by-P correlation matrix.  If U is an N-by-2
#'   matrix, RHO may also be a scalar correlation coefficient.
#'
#'   Y = COPULACDF(FAMILY,U,ALPHA) returns the cumulative probability of the
#'   bivariate Archimedean copula determined by FAMILY, with scalar parameter
#'   ALPHA, evaluated at the points in U.  FAMILY is 'Clayton', 'Frank', ort
#'   'Gumbel'.  U is an N-by-2 matrix of values in [0,1].
#'
#' @name  CopulaFamiliesCDF
#'
#' @examples  u = seq(0,1,0.1);
#'     U1=matrix(rep(u,length(u)),nrow=length(u),byrow = TRUE); U2=t(U1)
#'    F = CopulaFamiliesCDF('clayton',cbind(c(U1), c(U2)),1)
#'
#' @import mvtnorm
#'
#' @export
#' 
CopulaFamiliesCDF <- function(family,u,...){
  
  param= list(...)
  if( nargs()<3 ){stop("Wrong Number Of Inputs")}

  d=dim(u)[2]
  if( d<2 ){ stop("too few dimensions")}

  u[u<0] = 0
  u[u>1] = 1
  families = c("gaussian","t","clayton","frank","gumbel")
  if(!(family  %in% families)){ stop(paste(family,"is not valid as copula name"))}

  switch(family,
         "gaussian" = {
           alpha = param[[1]]

           if(d==2 & length(alpha==1)){
             Rho  = (2*exp(alpha)/(exp(alpha)+1)) - 1
             if(!(abs(Rho) < 1)){stop("Bad Scalar Correlation")}
             Rho = matrix(c(1,Rho, Rho,1),ncol=2,byrow=TRUE)
           } else if( dim(as.matrix(alpha))!=c(d,d)  | any(diag(alpha)!=Inf)){ stop("Bad Correlation Matrix")
           } else{
             Rho = ( 2*exp(alpha) / (exp(alpha)+1) ) - 1
             Rho[is.na(Rho)]=1
           }
           X=stats::qnorm(u)
           if(d==2){
             p = mvncdf(X , mean=rep(0,d) , Rho)
           } else {
             Rho = ( 2*exp(alpha) / (exp(alpha)+1) ) - 1
             Rho[is.na(Rho)]=1
             p = mvncdf(X , mean=rep(0,d) , Rho)
           }

         },

         "t" = {
           alpha  = param[[1]]
           alpha2 = param[[2]]
           if(d==2 & length(alpha==1)){
             Rho  = (2*exp(alpha)/(exp(alpha)+1)) - 1
             if(!( abs(Rho) < 1)){stop("Bad Scalar Correlation")}
             Rho = matrix(c(1,Rho, Rho,1),ncol=2,byrow=TRUE)
           }else if( dim(as.matrix(alpha))!=c(d,d)  | any(diag(alpha)!=Inf)){ stop("Bad Correlation Matrix")
           }else{
             Rho = ( 2*exp(alpha) / (exp(alpha)+1) ) - 1
             Rho[is.na(Rho)]=1
           }
           if( length(alpha2)!=1){stop("Bad Degrees Of Freedom")}
           p = mvtcdf(stats::qt(u,exp(alpha2)), Rho ,df = exp(alpha2))
         },

         "clayton" = {
           alpha = param[[1]]
           if(is.infinite(alpha)){ p = apply(u,1,prod)
           }else  p = (rowSums(u^(-exp(alpha))) - 1)^(-1/exp(alpha))
         },

         "frank" = {
           alpha = param[[1]]
           if(alpha==0){
             p = apply(u,1,prod)
           }else{
             p = -log((exp(-alpha) + (exp(-alpha*rowSums(u)) - rowSums(exp(-alpha*u)))) / expm1(-alpha)) / alpha
           }
         },

         "gumbel" = {
           alpha = param[[1]]
           if(is.infinite(alpha)){ p = apply(u,1,prod)
           }else  {
             v = matrix(Inf,ncol=dim(u)[2],nrow=dim(u)[1])
             for(j in 1:dim(v)[2]){
               v[u[,j]!=0,j]= -log(u[u[,j]!=0,j])
             }
             realmin = 2.2251e-308 ; realmax=  1.7977e+308
             # realmin and real max are Smallest and biggest positive normalized floating point number
             vmin = apply(u,1,min) ; vmax = apply(u,1,max)
             vmax[vmax==0] = realmin;  vmin[vmin==Inf] = realmax
             p = exp(-vmax*(1+(vmin/vmax)^(exp(alpha)+1))^ (1/(exp(alpha)+1)) )
           }
         }

  )
  return(p)
}





mvtcdf<- function(X,Rho,df){
  if(is.vector(X)){X = t(as.matrix(X))}

  p = rep(NA,dim(X)[1])
  for(i in 1:dim(X)[1]){
    p[i] = mvtnorm::pmvt(lower=rep(-Inf,dim(X)[2]),upper=X[i,],sigma=Rho,df=df)[1]
  }

  return(p)
}





mvncdf<-function(X,mean,Rho){
  if(is.vector(X)){X = t(as.matrix(X))}

  p = rep(NA,dim(X)[1])
  for(i in 1:dim(X)[1]){
    p[i] = mvtnorm::pmvnorm(lower=rep(-Inf,dim(X)[2]),upper=X[i,],mean=mean,sigma=Rho)[1]
  }

  return(p)
}


