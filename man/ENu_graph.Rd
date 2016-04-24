\name{ENu_graph}
\alias{ENu_graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots empirical expected number of upcrossings of level u with respect to P(Y<u) }
\description{
Plots empirical expected number of upcrossings of level u with respect to P(Y<u) }\usage{
ENu_graph(data, u, lty = 1, col = 1, add = FALSE, CI = FALSE, alpha = 0.05,
 N.s.data = NULL, xlab = "P(Y<u)", 
 ylab = "Intensity of upcrossings", ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.}
  \item{u}{
sequence of levels to be considered}
  \item{lty}{
type of line}
  \item{col}{
color of line}
  \item{add}{
if add=TRUE lines is added to current plot}
  \item{CI}{
if CI=TRUE a fluctuation interval at 1-alpha level of confidence is computed and plotted}
  \item{alpha}{
confidence level}
  \item{N.s.data}{
%%     ~~Describe \code{N.s.data} here~~
}
  \item{xlab}{
a title for the x axis}
  \item{ylab}{
a title for the y axis}
  \item{ylim}{
numeric vectors of length 2, giving the  y coordinates ranges.}
}

\value{list including
\item{u}{sequence of levels}
\item{F}{empirical cdf: P(data<u)}
\item{Nu}{number of upcrossings}
\item{CI.}{fluctuation interval}

}
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}



\seealso{
valid_all}
\examples{
data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
M = 2
order = 1
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH")
mod.hh= NULL
mod.hh$theta = theta.init
mod.hh$theta$A = matrix(c(0.40,0.88,-.09,-.13),2,2)
mod.hh$theta$A0 = matrix(c(6.75,1.08),2,1)
mod.hh$theta$sigma = matrix(c(1.76,3.40),2,1)
mod.hh$theta$prior = matrix(c(0.37,0.63),2,1)
mod.hh$theta$transmat = matrix(c(0.82,0.09,0.18,0.91),2,2)
#B.sim = 100*N.samples
#Y0 = array(data[1:2,sample(1:dim(data)[2],B.sim,replace=TRUE),],c(2,B.sim,1))
#Y.sim = simule.nh.MSAR(mod.hh$theta,Y0=Y0,T,N.samples=B.sim)
u = seq(min(data),max(data),by=.3)
gr.d = ENu_graph(data,u)
#gr = ENu_graph(Y.sim$Y,u,col=2,add=TRUE,CI = TRUE,N.s.data=dim(data)[2])
}

\keyword{ Mean upcrossings }
\keyword{ Validation }% __ONLY ONE__ keyword per line
