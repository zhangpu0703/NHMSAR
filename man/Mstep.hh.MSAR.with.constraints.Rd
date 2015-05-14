\name{Mstep.hh.MSAR.with.constraints}
\alias{Mstep.hh.MSAR.with.constraints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for fitting  homogeneous multivariate Markov switching auto-regressive models with constraints on VAR models.}
\description{
M step of the EM algorithm for fitting   homogeneous multivariate Markov switching auto-regressive models with constraints on VAR models, called in fit.MSAR. Maximum likelihood is used. Matrices A and sigma are diagonal by blocks.}
\usage{
Mstep.hh.MSAR.with.constraints(data, theta, FB, K, d.y)
}
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T x N.samples x d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{theta}{
model's parameter; object of class MSAR. See also init.theta.MSAR.
}
  \item{FB}{
Forward-Backward results, obtained by calling Estep.MSAR function}

  \item{K}{
number of sites. For instance, if one considers wind at k locations, K=k. Or more generally number of independent groups of components. }
  \item{d.y}{
dimension in each sites. For instance, if one considers only wind intensity than d.y = 1; but, if one considers cartesian components of wind, then d.y =2.}
}


\value{
	\item{A0}{intercepts}
	\item{A}{AR coefficients}
	\item{sigma}{variance of innovation}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}

}
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr
}


\seealso{
Mstep.hh.MSAR, fit.MSAR, Mstep.hh.SCAD.MSAR
}

