\name{emisprob.MSAR.VM}
\alias{emisprob.MSAR.VM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Emission probabilities for von Mises MSAR
	%%  ~~function to do ... ~~
}
\description{
Computes emission probabilities for von Mises MSAR models}
\usage{
emisprob.MSAR.VM(data, theta, covar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.}
  \item{theta}{
model's parameter; object of class MSAR. See also init.theta.MSAR.VM.}
  \item{covar}{
covariables for emission probabilities.}
}

\value{
prob : emission probabilities for each observation and each regime}
\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI. }
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
emisprob.MSAR}
