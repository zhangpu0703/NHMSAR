\name{Mstep.hh.SCAD.cw.MSAR}
\alias{Mstep.hh.SCAD.cw.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for fitting  homogeneous multivariate Markov switching auto-regressive models with SCAD penalization of parameters of the VAR(1) models.}
\description{
M step of the EM algorithm for fitting   homogeneous multivariate Markov switching auto-regressive models with penalization of parameters of the VAR(1) models, called in fit.MSAR. Penalization may be add to the autoregressive matrices of order 1 and to the precision matrices (inverse of variance of innovation). For the autoregressive matrices the ncvreg component wise procedure is used (see package ncvreg). For the precision matrices the graphcal lasso algortihm of glasso is used with the adaptative lasso of Zou.  }
\usage{
Mstep.hh.SCAD.cw.MSAR(data, theta, FB, lambda1=.1,lambda2=.1,penalty=,par=NULL)
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

  \item{lambda1}{
penalization constant for the precision matrices. It may be a scalar or a vector of length M (with M the number of regimes). If it is equal to0 no penalization is introduced for the precision matrices. }
  \item{lambda2}{
penalization constant for the autoregressive matrices. It may be a scalar or a vector of length M (with M the number of regimes).  If it is equal to0 no penalization is introduced for the atoregression matrices.}
\item{penalty}{choice of the penalty for the autoregressive matrices. Possible values are ridge, lasso or SCAD (default). }
\item{par}{allows to give an initial value to the precision matrices. }
}

\details{When LASSO penalty is chosen, the LARS algorithm is used. When SCAD is chosen, a Newton-Raphson algorithm is run with a quadratic approximation of the penalized likelihood. 
	For the precision matrices penalization, the package glasso is used. 
	
	Limit of this function: only works for VAR(1) models
		  }

\value{
	\item{A0}{intercepts}
	\item{A}{AR coefficients}
	\item{sigma}{variance of innovation}
	\item{sigma.inv}{inverse of variance of innovation}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}

}
\references{Breheny, P., & Huang, J. (2011). Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection. The annals of applied statistics, 5(1), 232.
	
	Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., et al. (2004). Least angle regression. The Annals of statistics, 32(2):407-499.

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. Journal of the American statistical Association, 96(456):1348-1360.

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9(3), 432-441.

}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr
}


\seealso{
Mstep.hh.MSAR, fit.MSAR, Mste.hh.SCAD.MSAR
}

