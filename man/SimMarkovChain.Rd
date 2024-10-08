% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimMarkovChain.R
\name{SimMarkovChain}
\alias{SimMarkovChain}
\title{Markov chain simulation}
\usage{
SimMarkovChain(Q, n, eta0)
}
\arguments{
\item{Q}{Transition probability matrix (d x d)}

\item{n}{number of simulated vectors}

\item{eta0}{variable eta}
}
\value{
\item{x}{Simulated Markov chain sequence}
}
\description{
Simulation of  n consecutive  values of a Markov chain with transition matrix Q, 
starting from a state eta0 or the uniform distribution on the set 1,..., r.
}
