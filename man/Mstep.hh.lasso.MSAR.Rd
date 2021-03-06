\name{Mstep.hh.lasso.MSAR}
\alias{Mstep.hh.lasso.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for fitting  homogeneous multivariate Markov switching auto-regressive models with penalization of parameters of the VAR(1) models.}
\description{
M step of the EM algorithm for fitting   homogeneous multivariate Markov switching auto-regressive models with penalization of parameters of the VAR(1) models, called in fit.MSAR. Penalized maximum likelihood is used. Penalization may be add to the autoregressive matrices of order 1 and to the precision matrices (inverse of variance of innovation).   }
\usage{
Mstep.hh.lasso.MSAR(data, theta, FB)
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

}

\details{The lars algorithm of pagkage lars is used.		  }

\value{
	\item{A0}{intercepts}
	\item{A}{AR coefficients}
	\item{sigma}{variance of innovation}
	\item{sigma.inv}{inverse of variance of innovation}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}

}
\references{Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., et al. (2004). Least angle regression. The Annals of statistics, 32(2):407-499.
}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr
}


\seealso{
Mstep.hh.MSAR, fit.MSAR
}

