\name{Mstep.hh.MSAR.VM}
\alias{Mstep.hh.MSAR.VM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for fitting  von Mises Markov switching auto-regressive models.}
\description{
M step of the EM algorithm for fitting   homogeneous Markov switching auto-regressive models, called in fit.MSAR.VM.}
\usage{
Mstep.hh.MSAR.VM(data, theta, FB, constr = 0)
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
Forward-Backward results, obtained by calling Estep.MSAR function
}
  \item{constr}{
constraints are added to the \eqn{\kappa} parameter (A preciser)}
}
\details{
The homogeneous MSAR model is labeled "HH" and it is written
 \deqn{
P(X_t|X_{t-1}=x_{t-1}) = Q_{x_{t-1},x_t}} with \eqn{X_t} the hidden univariate process defined on \eqn{\{1,\cdots,M \}}
 \deqn{
Y_t|X_t=x_t,y_{t-1},...,y_{t-p}} has a von Mises distribution with density \deqn{ p_2(y_t|x_t,y_{t-s}^{t-1}) = \frac {1}{b(x_t,y_{t-s}^{t-1})} \exp\left(\kappa_0^{(x_t)} \cos(y_t-\phi_0^{(x_t)})+
    \sum_{\ell=1}^s\kappa_\ell^{(x_t)} \cos(y_t-y_{t-\ell}-\phi_\ell^{(x)})\right)}
    which is equivalent to
\deqn{ p_2(y_t|x_t,y_{t-s}^{t-1}) = \frac {1}{b(x_t,y_{t-s}^{t-1})} \left|\exp\left([\gamma_0^{(x_t)} +
    \sum_{\ell=1}^s\gamma_\ell^{(x_t)} e^{iy_{t-\ell}}]e^{-iy_t}\right)\right|}
    
    \eqn{b(x_t,y_{t-s}^{t-1})} is a normalisation constant. 
    
    Both the real and the complex formulation are implemented. In practice, the complex version is used if the initial \eqn{\kappa} is complex. 
  }
\value{List containing
	\item{mu}{intercepts}
	\item{kappa}{von Mises AR coefficients}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}

}
\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{fit.MSAR.VM, Estep.MSAR.VM
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ M step }% __ONLY ONE__ keyword per line
