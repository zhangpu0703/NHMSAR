\name{forecast.prob.MSAR}
\alias{forecast.prob.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forecast probabilities for (non) homogeneous MSAR models}
\description{
Computes, for each time \eqn{t}, the conditional probabilities for MSAR models \eqn{P(Y_t|y_{1:(t-1)})} where  \eqn{Y} is the observed process and \eqn{y} the observed time series.}
\usage{
forecast.prob.MSAR(data, theta, yrange = NULL, covar.emis = NULL, covar.trans = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed time series, array of dimension T*N.samples*d 
}
  \item{theta}{
object of class MSAR including the model's parameter and description. See init.theta.MSAR for more details.}
  \item{yrange}{
values at which to compute the forecast probabilities
}
  \item{covar.emis}{
emission covariate if any.}
  \item{covar.trans}{
array of univariate or multivariate series of covariate to take into account in the transition probabilities.
The link function is defined in the initialisation step (running init.theta.MSAR.R).}
}

\value{A list containing
  \item{..$yrange }{abscissa for the forecast probabilities}
\item{..$prob }{forecast probabilities}
\item{Yhat}{forecasted value}
}
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
prediction.MSAR
}
\examples{
## Not run
#data(meteo.data)
#data = array(meteo.data$temperature,c(31,41,1)) 
#T = dim(data)[1]
#N.samples = dim(data)[2]
#d = dim(data)[3]
#M = 2
#theta.init = init.theta.MSAR(data,M=M,order=2,label="HH")
#res.hh.2 = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=200)
#FP = forecast.prob.MSAR(data,res.hh.2$theta)
#plot(data[,1,],typ="l")
#lines(FP$Yhat[,1],col="red")
#alpha = .1
#IC.emp = matrix(0,2,T)
#for (k in 1:length(data[,1,])) {
#	tmp = cumsum(FP$prob[,k,1])/sum(FP$prob[,k,1])
#	IC.emp[1,k] = FP$yrange[max(which(tmp<alpha/2))]
#	IC.emp[2,k] = FP$yrange[max(which(tmp<(1-alpha/2)))]
#}
#lines(IC.emp[1,],lty=2,col="red")
#lines(IC.emp[2,],lty=2,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Forecasting }
\keyword{ MSAR }% __ONLY ONE__ keyword per line
