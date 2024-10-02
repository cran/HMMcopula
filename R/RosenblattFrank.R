#'@title Rosenblatt transform for Frank copula
#'
#'@description Computation of the Rosenblatt transform for Frank's copula
#'
#'@param U   (n x d) matrix of pseudos-observations (normalized ranks)
#'@param theta parameter of the Frank copula
#'
#'@return \item{R}{Rosenblatt transform}
#'
#'@export
#'
RosenblattFrank <- function(U,theta){
  #  Rosenblatt transform for the Gaussian copula
  #  Input
  #        theta:  parameter  of the Frank;
  #         U : n x d matrix of pseudos-observations (normalized ranks).
  #
  #  By Bruno Remillard (Jul 27, 2015)


  CoefPolyFrank <- function(k){

    a = rep(0,k)
    for(j in 1:k){
      a[j] = polyFrank(k,j)
    }
    return(a)
  }

  polyFrank <- function(k,j){
    #  Computes polynomials associated with the Frank copula.

    if(k<=2){
      coeff = 1
    }else if(k>2){
      if(j==1){
        coeff = polyFrank(k-1,1)
      }else if( j==k){
        coeff = (k-1)*polyFrank(k-1,k-1)
      } else coeff = j*polyFrank(k-1,j)+ (j-1)*polyFrank(k-1,j-1)
    }
    return(coeff)
  }



  alpha = -log(theta) # copulafit matlab parameter
  n=dim(U)[1]; d=dim(U)[2]
  R = matrix(0,n,d)
  R[,1]=U[,1]
  x = -log(  (1-exp(-alpha*U))/(1-theta))
  v = 1 / (1-(1-theta)*exp(-t(apply(x,1,cumsum))))-1; #  exp(-alpha C)

  for(j in 2:d){
    l=j-1
    a = l/theta
    t1= v[,j]
    t2= v[,l]
    Coef =  CoefPolyFrank(l)
    p1 = rep(0,n)
    p2 = p1
    for(k in 1:l){
      p1 = Coef[k]* t1^k
      p2 = Coef[k]* t2^k
    }

    R[,j] =  p1 / p2
  }



  return(R)
}
