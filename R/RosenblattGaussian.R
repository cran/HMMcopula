#'@title Rosenblatt transform for Gaussian copula
#'
#' @description This function computes the Rosenblatt transform fot the Gaussian copula
#'
#'@param u   (n x d) matrix of pseudos-observations (normalized ranks)
#'@param rho (d x d) correlation matrix, or the correlation coefficient (if, d = 2)
#'
#'@return \item{R}{Rosenblatt transform}
#'
#'@export
#'
RosenblattGaussian <- function(u,rho){
  #  Rosenblatt transform for the Gaussian copula
  #  Input
  #        U: n x d matrix of pseudos-observations (normalized ranks);
  #        rho:  d x d correlation matrix.

  #
  #  By Bruno Remillard (Feb 10, 2010)
  n=dim(u)[1]; d=dim(u)[2]
  if (d == 2 & length(rho)==1 ){
    rho = matrix(c(1, rho, rho, 1),2)
  }
  R = matrix(0,n,d)
  R[,1]=u[,1]

  for(k in 2:d){

    R11 = rho[1:(k-1),1:(k-1)]
    R21 = rho[k,1:(k-1)]
    R11inv = solve(R11)
    B = R21%*%R11inv
    tu = stats::qnorm(u)
    x = tu[,1:(k-1)]
    y = tu[,k]
    Omega = 1-B%*%R11%*%t(B)
    mu = x%*%t(B)
    w = (y-mu) %*% solve(sqrt(Omega))
    R[,k] = stats::pnorm(w)

  }
  return(R)
}
