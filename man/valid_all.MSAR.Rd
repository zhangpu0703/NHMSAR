\name{valid_all.MSAR}
\alias{valid_all.MSAR}
\title{
Statistics plotting for validation of MSAR models
}
\description{
plots some functional statistics to help to valid  MSAR models: qqplot, covariance function, mean duration of sojourn over and under a threshold. For each of them the empirical statistic of the observed time series is plotted as well as the simulated one with \eqn{(1-\alpha)}-fluctuation intervals.}
\usage{
valid_all.MSAR(data,simu,title="",id=1,alpha=.05,spaghetti=TRUE,
mfrow=NULL,save=FALSE,output=FALSE,
root.filename=" ",path=NULL,col="red",width=4,height=4)
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
\item{spaghetti}{statistics of every simulation batch are plotted instead of fluctuation intervals. A batch is a simulation block of the same size as the observations. Default spaghetti=TRUE}
 \item{mfrow}{if NULL, each plot is done in a new window}
 \item{save}{if save=TRUE plots are saved into .eps files}
  \item{root.filename}{
root file name for saving plots}
  \item{path}{
path of folder where to save the files}
  \item{output}{
if TRUE some statistics are returned.  }
\item{col}{color of the lines for simulated data, default is red}
\item{width}{width of the figure when is it save by dev.copy2eps}
\item{height}{height of the figure when is it save by dev.copy2eps}
}

\value{Returns plots and
	\item{qqp}{statistics of marginal distributions}
	\item{C}{statistics of correlation functions}
	\item{ENu.data}{statistics of  intensity of up crossings of the data}
	\item{ENu.simu}{statistics of  intensity of up crossings of the simulations}
	\item{MDO}{statistics of  mean duration over a level}
	\item{MDU}{statistics of  mean duration under a level}
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
valid_all.MSAR(data,Y.sim$Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MSAR model validation }
