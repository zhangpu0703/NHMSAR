\name{Mstep.nh.MSAR}
\alias{Mstep.nh.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm.}
\description{
M step of the EM algorithm for fitting Markov switching auto-regressive models with non homogeneous transitions.}
\usage{
Mstep.nh.MSAR(data, theta, FB, covar = NULL, method = method,ARfix=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{theta}{
model's parameter; object of class MSAR. See also init.theta.MSAR.
}
  \item{FB}{
Forward-Backward results, obtained by calling Estep.MSAR function}

  \item{covar}{
transitions covariates}
  \item{method}{
permits to choice the optimization algorithm. default is "ucminf", other possible choices are "BFGS" or "L-BFGS-B"}
 \item{ARfix}{
if TRUE the AR parameters are not estimated, they stay fixed at their initial value. }
}

\value{
List containing
	\item{..$A0}{intercepts}
	\item{..$A}{AR coefficients}
	\item{..$sigma}{variance of innovation}
	\item{..$prior}{prior probabilities}
	\item{..$transmat}{transition matrix}
	\item{..$par.trans}{transitions parameters}

}

\references{
Ailliot P., Monbet V., (2012), Markov switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR, init.theta.MSAR, Mstep.hh.MSAR}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EM algorithm }
\keyword{ M step }
