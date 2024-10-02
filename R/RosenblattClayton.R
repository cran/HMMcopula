#'@title Rosenblatt transform for Clayton copula
#'
#'@description Computation of the Rosenblatt transform fot Clayton's copula
#'
#'@param u   (n x d) matrix of pseudos-observations (normalized ranks)
#'@param theta parameter of the Clayton copula
#'
#'@return \item{R}{Rosenblatt transform}
#'
#'@export
#'
RosenblattClayton <- function(u,theta){
  #  Rosenblatt transform for the Gaussian copula
  #  Input
  #        U: n x d matrix of pseudos-observations (normalized ranks);
  #        theta: parameter of the Clayton;

  #
  #  By Bruno Remillard (Feb 11, 2010)
  n=dim(u)[1]; d=dim(u)[2]
  R = matrix(0,n,d)
  R[,1]=u[,1]
  x = u^(-theta)-1
  u = 1+ t(apply(x,1,cumsum))

  for(j in 2:d){

    a = -1/theta-(j-1)
    z1 = (u[,j])^a
    z2 = (u[,j-1])^a
    R[,j] = z1/z2

  }
  return(R)
}
