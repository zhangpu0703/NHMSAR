\name{MeanDurOver}
\alias{MeanDurOver}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Duration of sojourn over a treshold
}
\description{
Plot the mean duration of sojourn over thresholds for an observed time series and a simulated one with respect to the empirical cumulative distribution function. Fluctuation intervals are plotted too.}
\usage{
MeanDurOver(data, data.sim, u, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*1 }
  \item{data.sim}{
simulated time series, array of dimension T*N.sim*1.
N.sim have to be  K*N.samples with K large enough (for instance, K=100) }
  \item{u}{vector of 
thresholds}
  \item{alpha}{1-confidence level for fluctuation intervals. Default = 0.05}}

\value{
Returns a plot and a list including
..$F : empirical cdf of data for levels u
..$mdo.data : mean duration over levels u for data
..$F.sim : empirical cdf of simulations for levels u
..$mdo.sim : mean duration over levels u for simulations
..$CI : confidence intervals of mean duration over levels u for simulations
..$mod.sim.all : mean duration over levels u for all simulations
}
%\references{
%% ~put references to the literature/web site here }
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{valid_all}
\examples{
data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
M = 2
order = 2
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH")
mod.hh= NULL
mod.hh$theta = theta.init
mod.hh$theta$A = matrix(c(0.40,0.88,-.09,-.13),2,2)
mod.hh$theta$A0 = matrix(c(6.75,1.08),2,1)
mod.hh$theta$sigma = matrix(c(1.76,3.40),2,1)
mod.hh$theta$prior = matrix(c(0.37,0.63),2,1)
mod.hh$theta$transmat = matrix(c(0.82,0.09,0.18,0.91),2,2)
B.sim = 50*N.samples
Y0 = array(data[1:2,sample(1:dim(data)[2],B.sim,replace=TRUE),],c(2,B.sim,1))
Y.sim = simule.nh.MSAR(mod.hh$theta,Y0=Y0,T,N.samples=B.sim)
u = seq(min(data),max(data),length.out=30)
MDO = MeanDurOver(data,Y.sim$Y,u)
}

\keyword{ Mean Duration of Sojourn }
\keyword{ Threshold excess }% __ONLY ONE__ keyword per line
