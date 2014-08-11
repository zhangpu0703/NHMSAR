\name{cross.cor.MSAR}
\alias{cross.cor.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
empirical cross-correlation for multivariate MSAR time series }
\description{
cross.cor.MSAR computes the cross-correlation between two components. The cross-corelation can be estimted for the whole time series or regime by regime. }
\usage{
cross.cor.MSAR(data, X=NULL, nc1 = 1, nc2 = 2, lag = 10, regime = 0, CI = FALSE, 
  Bsim = 0, N.samples = 1, add = FALSE, col = 1, names = NULL, alpha = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*d }
  \item{X}{
time series of regimes associated to data}
  \item{nc1}{
first component to be considered}
  \item{nc2}{
second component to be considered}
  \item{lag}{
maximum lag (default=10). The cross-correlation is estimated for lags -lag:lag.
}
  \item{regime}{
has to be an integer between 0 and M, with M the number of regimes. If regime=0, the cross correlaiton is computed for the whole time series. If regime=m>0, the corss corelation is computed considereing only the sub-sequences in regime m.  }
  \item{CI}{
If CI=TRUE fluctuation intervals are computed, default is FALSE}
  \item{Bsim}{
useful for computation of confidence intervals. When observed and simulated data are compared, one expects that the number of simulated time series is Bsim*N.samples }
  \item{N.samples}{
useful for computation of confidence intervals. N.sample describes the number of independant time series in the observed (or reference) data }
  \item{add}{
if add=TRUE the empirical cross-correlation is added to the current plot.}
  \item{col}{
color of the line}
  \item{names}{
list with the names of components of data}
  \item{alpha}{
level for the computation of the fluctuation intervals. default=0.1}
}
\details{The cross-correlation functions are computed from one or several independent realizations of the same length.
%%  ~~ If necessary, more details than the description above ~~
}
\value{returns a list including:
	\item{..$ccf}{empirical cross-correlation}
	\item{..$lag}{abscissa for the cross-correlation}
	\item{..$CI}{fluctuation intervals}


}
\references{
Bessac, J., Ailliot, P., & Monbet, V. (2013). Gaussian linear state-space model for wind fields in the North-East Atlantic. arXiv preprint arXiv:1312.5530.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


\seealso{
cor.MSAR, cor, valid_all}
\examples{
data(Wind)
T = dim(U)[1]
c = cross.cor.MSAR(U,nc1=1,nc2=18,names=1:18)
## Not run
#Y = U[,,c(1,18)]
#theta.init=init.theta.MSAR(Y,M=2,order=2,label="HH")
#res.hh = fit.MSAR(Y,theta.init,verbose=TRUE,MaxIter=200)
#Bsim = 20
#N.samples = dim(U)[2]
#Ksim = Bsim*N.samples
#Y0 = Y0
#Y.sim = simule.nh.MSAR(res.hh$theta,Y0 = Y0,T,N.samples = Ksim) 
#c.sim = cross.cor.MSAR(Y.sim$Y,nc1=1,nc2=2,names=c(1,18),
#   CI=TRUE,Bsim=Bsim,N.samples=N.samples,add=TRUE,col="red")
}

\keyword{ Cross-correlation function }
\keyword{ Validation }% __ONLY ONE__ keyword per line
