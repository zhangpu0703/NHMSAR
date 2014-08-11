\name{Estep.MSAR.VM}
\alias{Estep.MSAR.VM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estep of the EM algorithm for fitting von Mises (non) homogeneous Markov switching auto-regressive models.}
\description{
Forward-backward algorithm called in fit.MSAR.}
\usage{
Estep.MSAR.VM(data, theta, smth = FALSE, verbose = FALSE, 
   covar.emis = NULL, covar.trans = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{theta}{
model's parameter; object of class MSAR. See also init.theta.MSAR.}
  \item{smth}{
If smth=FALSE, only the forward step is computed for forecasting probabilities. If smth=TRUE, the smoothing probabilities are computed too. }
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
   \item{covar.emis}{
covariables for emission probabilities.}
  \item{covar.trans}{
covariables for transition probabilities}
}
\value{
list including
\item{loglik}{log likelihood}
\item{probS}{smoothing probabilities: \eqn{P(S_t=s|y_0,\cdots,y_T)}}
\item{probSS}{one step smoothing probabilities: \eqn{P(S_t=s,S_{t+1}|y_0,\cdots,y_T)}}

}
\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

\seealso{
fit.MSAR.VM, Mstep.hh.MSAR.VM,Estep.MSAR}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ E step }% __ONLY ONE__ keyword per line
\keyword{ Forward-backward }
