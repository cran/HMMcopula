% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimHMMCop.R
\name{SimHMMCop}
\alias{SimHMMCop}
\title{Simulation of bivariate Markov regime switching copula model}
\usage{
SimHMMCop(Q, family, KendallTau, n, DoF)
}
\arguments{
\item{Q}{Transition probality matrix (d x d);}

\item{family}{'gaussian' , 't' , 'clayton' , 'frank' , 'gumbel'}

\item{KendallTau}{Kendall's rank correlation}

\item{n}{number of simulated vectors}

\item{DoF}{degree of freedom only for the Student copula}
}
\value{
\item{SimData}{Simulated Data}

\item{MC}{Markov chain regimes}

\item{alpha}{parameters alpha}
}
\description{
Simulaion of values from a bivariate Markov regime switching copula model
}
\examples{
Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) ; kendallTau <- c(0.3 ,0.7) ;
simulations <- SimHMMCop(Q, 'gumbel', kendallTau, 300)


}
