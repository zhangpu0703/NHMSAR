\name{MeanDurUnder}
\alias{MeanDurUnder}
\title{
Mean Duration of sojourn under a treshold
}
\description{
Plot the mean duration of sojourn under thresholds for an observed time series and a simulated one with respect to teh empirical cumulative distribution function (cdf).  Confidence intervals are plotted too.}

\usage{
MeanDurUnder(data, data.sim, u, alpha = 0.05,col="red",plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*1 }
  \item{data.sim}{
simulated time series, array of dimension T*N.sim*1.
N.sim have to be  K*N.samples with K large enough (for instance, K=100) }
  \item{u}{
vector of thresholds}
  \item{alpha}{
1-confidence level for confidence intervals. Default = 0.05}
\item{col}{color of the lines for simulated data, default is red}
\item{plot}{statistic are plotted if TRUE (default)}
}

\value{
Returns a plot and a list including
..$F : empirical cdf of data for levels u
..$mdu.data : mean duration under levels u for data
..$F.sim : empirical cdf of simulations for levels u
..$mdu.sim : mean duration under levels u for simulations
..$CI : confidence intervals of mean duration under levels u for simulations
}
%\references{
%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

\seealso{valid_all.MSAR, MeanDurOver}
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
B.sim = 20*N.samples
Y0 = array(data[1:2,sample(1:dim(data)[2],B.sim,replace=TRUE),],c(2,B.sim,1))
Y.sim = simule.nh.MSAR(mod.hh$theta,Y0=Y0,T,N.samples=B.sim)
u = seq(min(data),max(data),length.out=30)
MeanDurUnder(data,Y.sim$Y,u)
}

\keyword{ Mean Duration of Sojourn }
\keyword{ Threshold excess }
