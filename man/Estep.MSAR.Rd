\name{Estep.MSAR}
\alias{Estep.MSAR}
\title{
Estep of the EM algorithm for fitting (non) homogeneous Markov switching auto-regressive models. }
\description{
Forward-backward algorithm called in fit.MSAR. }
\usage{
Estep.MSAR(data, theta, smth = FALSE, 
           verbose = FALSE, 
           covar.emis = covar.emis, covar.trans = covar.trans)
}
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{theta}{model's parameter; object of class MSAR. See also init.theta.MSAR.
.}
  \item{smth}{
If smth=FALSE, only the forward step is computed for forecasting probabilities. 
If smth=TRUE, the smoothing probabilities are computed too. }
  \item{verbose}{
if verbose=TRUE some results are printed at each iteration. }
  \item{covar.emis}{
covariables for emission probabilities.}
  \item{covar.trans}{
covariables for transition probabilities.}
}

\value{
A list including
\item{loglik}{log likelihood}
\item{probS}{smoothing probabilities: \eqn{P(S_t=s|y_0,\cdots,y_T)}}
\item{probSS}{one step smoothing probabilities: \eqn{P(S_t=s,S_{t+1}|y_0,\cdots,y_T)}}

}
\references{
Ailliot P., Monbet V., (2012), Markov switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}




\seealso{fit.MSAR, Mstep.hh.MSAR}
\examples{#see fit.MSAR}

\keyword{ EM algorithm }
\keyword{ E step }
\keyword{ Forward-backward }