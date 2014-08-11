\name{valid_all}
\alias{valid_all}
\title{
Statistics plotting for validation of MSAR models
}
\description{
plots some functional statistics to help to valid  MSAR models: qqplot, covariance function, mean duration of sojourn over and under a threshold. For each of them the empirical statistic of the observed time series is plotted as well as the simulated one with \eqn{(1-\alpha)}-fluctuation intervals.}
\usage{
valid_all(data, simu, root.filename = " ", path = NULL, title = "", 
  id = 1, alpha = 0.05, save = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*d }
  \item{simu}{
  simulated time series, array of dimension T*N.sim*d.
N.sim have to be  K*N.samples with K large enough (for instance, K=100) }
   \item{title}{
title of plots}
  \item{id}{
component to be considered when the data is multivariate (d>1). Default d=1.}
  \item{alpha}{
level for the \eqn{(1-\alpha)}-fluctuation intervals}
 \item{save}{
if save=TRUE plots are saved into .eps files}
  \item{root.filename}{
root file name for saving plots}
  \item{path}{
path of folder where to save the files}
}

\value{Returns plots 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}




\examples{
	data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
plot(data[,k,1],typ="l",xlab=("time (days)"),ylab=("temperature (degrees C)"))
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
# Fit Homogeneous MS-AR models
M = 2
order = 1
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH") 
mod.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=10)
# Simulation
yT = 31
Bsim = 10
Ksim = Bsim*N.samples
Y0 = array(data[1:2,sample(1:dim(data)[2],Ksim,replace=T),],c(2,Ksim,1))
Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = Ksim) 
valid_all(data,Y.sim$Y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MSAR model validation }
