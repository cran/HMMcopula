#'@title COPULAPDF Probability density function for a copula.
#'
#'COPULAPDF  probability density function for a copula with linear correlation parameters RHO and
#'
#'@param family  copula familly= "gaussian" , "t" , "clayton" , "frank" , "gumbel"
#'@param u  is an N-by-P matrix of values in [0,1], representing N points in the P-dimensional unit hypercube
#'@param ...  additionnal parameter like RHO  a P-by-P correlation matrix.
#'
#'@return
#'Y = COPULAPDF('Gaussian',U,RHO) returns the probability density of the
#'   Gaussian copula with linear correlation parameters RHO, evaluated at the
#'   points in U. U is an N-by-P matrix of values in [0,1], representing N
#'   points in the P-dimensional unit hypercube.  RHO is a P-by-P correlation
#'   matrix.  If U is an N-by-2 matrix, RHO may also be a scalar correlation
#'   coefficient.
#'
#'   Y = COPULAPDF('t',U,RHO,NU) returns the probability density of the t
#'   copula with linear correlation parameters RHO and degrees of freedom
#'   parameter NU, evaluated at the points in U.  U is an N-by-P matrix of
#'   values in [0,1]. RHO is a P-by-P correlation matrix.  If U is an N-by-2
#'   matrix, RHO may also be a scalar correlation coefficient.
#'
#'   Y = COPULAPDF(FAMILY,U,ALPHA) returns the probability density of the
#'   bivariate Archimedean copula determined by FAMILY, with scalar parameter
#'   ALPHA, evaluated at the points in U.  FAMILY is 'Clayton', 'Frank', or
#'   'Gumbel'.  U is an N-by-2 matrix of values in [0,1].
#'
#'@examples  u = seq(0,1,0.1);
#'     U1=matrix(rep(u,length(u)),nrow=length(u),byrow = TRUE); U2=t(U1)
#'    F = copulaFamiliesPDF('clayton',cbind(c(U1), c(U2)),1)
#'
#'@import mvtnorm matrixcalc foreach doParallel copula
#'
#'
#'@export
#'@keywords internal
copulaFamiliesPDF<-function(family,u,...){
             param= list(...)
             if( nargs()<3 ){stop("Wrong Number Of Inputs")}
             n = dim(u)[1]
             d = dim(u)[2]

             if( d<2 ){ stop(" too few dimensions")}
             outofRange = apply(u<=0 | 1<=u,1,any)
             if(any(outofRange)){
               # Replace entire rows that include out-of-range values with 0.5, but
               # don't overwrite NaN. This will help protect against NaN caused by
               # data close to the edge of the space. Any NaN found in the output
               # should then be the result of NaN input or bad parameter values.
               for(i in 1:dim(u)[1]){
                 if( all(!is.na(u[i,])) & any(u[i,]<=0 | 1<=u[i,]) ){
                   u[i,]=0.5
                 }
               }
             }

             families = c("gaussian","t","clayton","frank","gumbel")
             if(!(family  %in% families)){ stop(paste(family,"is not valid as copula name"))}

             switch(family,
                    "gaussian" = {
                      alpha = param[[1]]
                      alpha = pmin(alpha,log(19999))
                      alpha = max(alpha,-log(19999))
                      if(d==2 & length(alpha==1)){
                        Rho  = (2*exp(alpha)/(exp(alpha)+1)) - 1
                        if(!(-1 < Rho & Rho < 1)){stop("Bad Scalar Correlation")}
                        Rho = matrix(c(1,Rho, Rho,1),ncol=2,byrow=TRUE)
                      } else if( dim(as.matrix(alpha))!=c(d,d)  | any(diag(alpha)!=Inf)){ stop("Bad Correlation Matrix")
                      } else{
                        Rho = ( 2*exp(alpha) / (exp(alpha)+1) ) - 1
                        Rho[is.na(Rho)]=1
                      }
                      R<-try(chol(Rho),silent = TRUE)
                      res <- try(matrixcalc:: is.positive.definite(Rho),silent = TRUE)
                      if(class(res) == "try-error"){ stop("Singular Correlation Matrix for Rho")}
                      x = stats::qnorm(u)
                      logSqrtDetRho = sum(log(diag(R)))
                      z = x %*% solve(R)
                      y =  exp(-0.5 * rowSums(z^2 - x^2) - logSqrtDetRho)
                    },

                    "t" = {
                      alpha  = param[[1]]
                      alpha = pmin(alpha,log(19999))
                      alpha = max(alpha,-log(19999))
                      alpha2 = param[[2]]
                      if(d==2 & length(alpha==1)){
                        Rho  = (2*exp(alpha)/(exp(alpha)+1)) - 1
                        if(!(-1 < Rho & Rho < 1)){stop("Bad Scalar Correlation")}
                        Rho = matrix(c(1,Rho, Rho,1),ncol=2,byrow=TRUE)
                      }else if( dim(as.matrix(alpha))!=c(d,d)  | any(diag(alpha)!=Inf)){ stop("Bad Correlation Matrix")
                      }else{
                        Rho = ( 2*exp(alpha) / (exp(alpha)+1) ) - 1
                        Rho[is.na(Rho)]=1
                      }
                      R<-try(chol(Rho),silent = TRUE)
                      res <- try(matrixcalc:: is.positive.definite(Rho),silent = TRUE)
                      if(class(res) == "try-error"){ stop("Singular Correlation Matrix for Rho")}


                      if(length(alpha2)<1){ stop("Missing Degree of freedom")}
                      if( length(alpha2)!=1){stop("Bad Degrees Of Freedom")}
                      tq = stats::qt(u,exp(alpha2))
                      z =  tq %*% solve(R)
                      logSqrtDetRho = sum(log(diag(R)))
                      const = lgamma((exp(alpha2)+d)/2) + (d-1)*lgamma(exp(alpha2)/2) - d*lgamma((exp(alpha2)+1)/2) - logSqrtDetRho
                      numer = -((exp(alpha2)+d)/2) * log(1 + rowSums(z^2)/exp(alpha2))
                      denom = rowSums(-((exp(alpha2)+1)/2) * log(1 + (tq^2)/exp(alpha2)))
                      y = exp(const + numer - denom)
                      },

                    "clayton" = {
                      # a.k.a. Cook-Johnson
                      # alpha = ln(du aplha matlab) = ln(du theta Bruno)
                      # alpha de matlab dans (0, +Inf)
                      # matlab ne prend pas en compte un (theta Bruno) dans (-1, 0)
                      # de ce fait alpha (-Inf , +Inf)
                      alpha = param[[1]]
                      if(is.infinite(alpha)){ y=rep(1,n)
                      }else {
                        logC = (-1/exp(alpha))* log(rowSums(u^(-exp(alpha))) - 1)
                        y = (exp(alpha)+1) * exp((2*exp(alpha)+1)*logC - rowSums((exp(alpha)+1)*log(u)))
                      }

                    },

                    "frank" = {

                      # alpha = alpha de matlab = -ln(theta de Bruno)
                      # alpha de matlab dans (-Inf, +Inf)
                      # de ce fait alpha (-Inf, +Inf)
                      alpha = param[[1]]
                      if(alpha==0){
                        y=rep(1,n)
                      }else{
                        y = alpha*(1-exp(-alpha))/ (cosh(alpha * t(apply(u,1,diff))/2)*2 - exp(alpha*(rowSums(u)-2)/2) - exp(-alpha*rowSums(u)/2))^2
                      }
                    },

                    "gumbel" = {
                      # a.k.a. Gumbel-Hougaard
                      # C(u1,u2) = exp(-((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))
                      # alpha de matlab = 1/(theta de Bruno);
                      # alpha matlab (1 , +Inf)
                      # alpha = ln(alpha de matlab - 1) = ln(1/theta - 1);
                      # de ce fait alpha (-Inf , +Inf)
                      alpha = param[[1]]
                      if(is.infinite(alpha)){
                        y=rep(1,n)
                      }else  {
                        v = -log(u)
                        vmin = apply(v,1,min) ; vmax = apply(v,1,max)
                        nlogC = vmax*(1+(vmin/vmax)^(exp(alpha)+1))^(1/(exp(alpha)+1))
                        y = (exp(alpha) + nlogC) * exp(-nlogC + rowSums((exp(alpha))*log(v) + v) + (1-2*(exp(alpha)+1))*log(nlogC))
                      }
                      if( d<2 ){ stop(" too few dimensions")}
                      outofRange = apply(u<=0 | 1<=u,1,any)
                    }

             )

             if(any(outofRange )){
               y[!is.na(y) & outofRange]=0
             }
             return(y)
           }



