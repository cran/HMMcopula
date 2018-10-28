#'Sample Kendall's tau Estimation
#'
#'@description This function estimates the sample Kendall's tau of a bivariate data matrix
#'@param X (n x 2) matrix
#'
#'@return \item{KendallTau}{estimated sample Kendall's tau of the data}
#'
#'@export
#'
EstKendallTau<-function(X){

  Nb_paire_concordante=0
  Nb_paire_dis_concordante=0

  for(ii in 1:(dim(X)[1]-1)){
    for(jj in (ii+1):dim(X)[1]){
      if(  ( (X[ii,1]-X[jj,1]) * (X[ii,2]-X[jj,2]) ) > 0  ) {
        Nb_paire_concordante=Nb_paire_concordante+1
      }

      if(  ( (X[ii,1]-X[jj,1]) * (X[ii,2]-X[jj,2]) ) < 0  ) {
        Nb_paire_dis_concordante= Nb_paire_dis_concordante+1
      }
    }
  }
  KendallTau = (Nb_paire_concordante - Nb_paire_dis_concordante)/(Nb_paire_concordante + Nb_paire_dis_concordante)
  return(KendallTau)
}
