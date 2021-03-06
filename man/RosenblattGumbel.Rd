% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RosenblattGumbel.R
\name{RosenblattGumbel}
\alias{RosenblattGumbel}
\title{Rosenblatt transform for Gumbel copula}
\usage{
RosenblattGumbel(U, theta)
}
\arguments{
\item{U}{(n x d) matrix of pseudos-observations (normalized ranks)}

\item{theta}{parameter of the Gumbel copula}
}
\value{
\item{R}{Rosenblatt transform}
}
\description{
This function computes the Rosenblatt transform fot the Gumbel copula
}
