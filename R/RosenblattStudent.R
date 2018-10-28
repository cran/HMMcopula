#'@title Rosenblatt transform for Student copula
#'
#' @description This function computes the Rosenblatt transform fot the Student copula
#'
#'@param u   (n x d) matrix of pseudos-observations (normalized ranks)
#'@param rho (d x d) correlation matrix
#'@param nu  degrees of freedom
#'
#'@return \item{R}{Rosenblatt transform}
#'
#'@export
#'
RosenblattStudent <- function(u,rho,nu){
  #  Rosenblatt transform for the Gaussian copula
  #  Input
  #        U: n x d matrix of pseudos-observations (normalized ranks);
  #        rho:  d x d correlation matrix.
  #        nu: degrees of freedom.
  #
  #  By Bruno Remillard (Feb 05, 2010)
  n=dim(u)[1]; d=dim(u)[2]

  if(d == 2 & length(rho)==1 ){
    rho = matrix(c(1, rho, rho, 1),2) }

  R = matrix(0,n,d)
  R[,1]=u[,1]

  for(k in 2:d){

    R11 = rho[1:(k-1),1:(k-1)]
    R21 = rho[k,1:(k-1)]
    R11inv = solve(R11)
    B = R21%*%R11inv
    tu = stats::qt(u, df = nu)
    x = tu[,1:(k-1)]
    y = tu[,k]
    Omega = 1-B%*%R11%*%t(B)
    mu = x%*%t(B)
    z = diag( nu + x%*%R11inv%*%t(x))/(nu+k-1)
    w = c(y-mu)/t(sqrt(Omega %*% z))
    R[,k] = stats::pt(w, df= nu+k-1)

  }
  return(R)
}
