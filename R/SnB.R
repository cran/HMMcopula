#'@title Cramer-von Mises statistic SnB for GOF based on the Rosenblatt transform
#'
#'@description This function computes the Cramer-von Mises statistic SnB for GOF based on the Rosenblatt transform
#'
#'@param E   (n x d) matrix of pseudos-observations (normalized ranks)
#'
#'@return \item{Sn}{Cramer-von Mises statistic}
#'
#'@export
#'
SnB <- function(E){
  #  Cramer-von Mises statistic SnB for GOF based on the Rosenblatt transform
  #  Ref: Genest, Remillard & Beaudoin 2009
  #
  #   Input
  #        E: n x d matrix of pseudo-observations.
  #
  #   Output
  #        Sn: Cramer-von Mises statistic SnB.

  n = dim(E)[1]; d = dim(E)[2]
  Dn= rep(0,n)
  S1= n/3^d
  G0 = 1-E*E
  E0 = 1-E
  S2 = sum(apply(G0,1,prod))/2^(d-1)

  for(i in 1:n){
    G0= matrix(rep(E0[i,],n), nrow=n, byrow = T)
    Dn[i]= mean( apply(pmin(G0,E0),1,prod))
  }
  Sn = sum(Dn)-S2+S1

  return(Sn)
}
