\name{Mstep.nn.MSAR}
\alias{Mstep.nn.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm.}
\description{
M step of the EM algorithm for fitting Markov switching auto-regressive models with non homogeneous emissions and non homogeneous transitions.}
\usage{
Mstep.nn.MSAR(data, theta, FB, 
   covar.trans = covar.trans, covar.emis = covar.emis, method = NULL)
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

  \item{covar.trans}{
transitions covariates}
  \item{covar.emis}{
emissions covariates (the covariates act on the intercepts)}
  \item{method}{
permits to choice the optimization algorithm. default is "ucminf", other possible choices are "BFGS" or "L-BFGS-B}
}

\value{
	\item{A0}{intercepts}
	\item{A}{AR coefficients}
	\item{sigma}{variance of innovation}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}
	\item{par_emis}{emission parameters}
	\item{par.trans}{transitions parameters}
}
\references{
Ailliot P., Monbet V., (2012), Markov switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Mstep.hh.MSAR}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ M step }% __ONLY ONE__ keyword per line
