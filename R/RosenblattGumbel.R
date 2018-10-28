#'@title Rosenblatt transform for Gumbel copula
#'
#'@description This function computes the Rosenblatt transform fot the Gumbel copula
#'
#'@param U   (n x d) matrix of pseudos-observations (normalized ranks)
#'@param theta parameter of the Gumbel copula
#'
#'@return \item{R}{Rosenblatt transform}
#'
#'@export
#'
RosenblattGumbel <- function(U,theta){
  #  Rosenblatt transform for the Gaussian copula
  #  Input
  #        theta: parameter of the Gumbel; Matlab parameter is 1/theta
  #         U : n x d matrix of pseudos-observations (normalized ranks).
  #
  #  By Bruno Remillard (Jul 27, 2015)
  n=dim(U)[1]; d=dim(U)[2]
  R = matrix(0,n,d)
  R[,1]=U[,1]
  x =  (-log(U))^(1/theta)
  u = 1+ t(apply(x,1,cumsum))
  v = t(apply(x,1,cumsum))^theta; # -log(C)

  CoefPolyGumbel <- function(k,theta){
    a = rep(0,k)
    for(j in 1:k){
      a[j] = polyGumbel(k,theta,j)
    }

    return(a)
  }

  polyGumbel <- function(k,theta,j){
    # Computes polynomials associated with the Gumbel copula of parameter theta
    # in [0,1].

    if(k==1){
      coeff = theta
    }else if(k==2){
      if(j==1){
        coeff = theta*(1-theta)
      }else if( j==2){
        coeff = theta^2
      }
    }else if(k>2){
      if(j==1){
        coeff = (k-1-theta)*polyGumbel(k-1,theta,1)
      }else if(j==k){
        coeff = theta*polyGumbel(k-1,theta,k-1)
      }else coeff = (k-1-j*theta)*polyGumbel(k-1,theta,j)+theta*polyGumbel(k-1,theta,j-1)

    }

    return(coeff)

  }



  for(j in 2:d){
    l=j-1
    a = l/theta
    t1= v[,j]
    t2= v[,l]
    Coef =  CoefPolyGumbel(l,theta)
    p1 = rep(0,n)
    p2 = p1
    for(k in 1:l){
      p1 = Coef[k]* t1^k
      p2 = Coef[k]* t2^k
    }

    R[,j] = exp(t2-t1) * (t2 / t1)^a * (p1 / p2)


  }



  return(R)
}


