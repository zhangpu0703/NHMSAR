\name{Mstep.hn.MSAR}
\alias{Mstep.hn.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for fitting Markov switching auto-regressive models with non homogeneous emissions.}
\description{
The M step contains two parts. One for the estimation of the parameters of the hidden Markov chain and the other for the parameters of the auto-regressive models. A numerical algortihm is used for the emission parameters. }
\usage{
Mstep.hn.MSAR(data, theta, FB, covar = NULL, verbose = FALSE)
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
emissions covariates (the covariables act on the intercepts)}
 \item{verbose}{if verbose is TRUE some iterations of the numerical optimisation are print on the console. }
}
\details{
The default numerical optimization method is ucminf (see ucminf). %%  ~~ If necessary, more details than the description above ~~
}
\value{
List containing
	\item{..$A0}{intercepts}
	\item{..$A}{AR coefficients}
	\item{..$sigma}{variance of innovation}
	\item{..$prior}{prior probabilities}
	\item{..$transmat}{transition matrix}
	\item{..$par_emis}{emission parameters}

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
\keyword{ M step }% __ONLY ONE__ keyword per line
