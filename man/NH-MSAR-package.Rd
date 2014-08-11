\name{NH-MSAR-package}
\alias{NH-MSAR-package}
\alias{NH-MSAR}
\docType{package}
\title{
 (Non) Homogeneous Markov switching autoregressive model
}
\description{
NH-MSAR-package is a set of functions to fit, simulate and validate (non) homogeneous Markov Switching Autoregressive models with Gaussian or von Mises innovations.}
\details{
\tabular{ll}{
Package: \tab NH-MSAR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-11\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Val\'e'rie Monbet, valerie.monbet@univ-rennes1.fr }
\references{
Hamilton J.D. (1989). A New Approach to the Economic Analysis of Nonstionary Time Series and
the Business Cycle. Econometrica 57: 357-384.
Ailliot P., Monbet V., (2012), Markov-switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101.
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.
}
\keyword{ package }

\examples{
	# Fit Homogeneous MS-AR models - univariate time series
data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
M = 2
order = 2
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH") 
mod.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=20)
regimes.plot.MSAR(mod.hh,data,ylab="temperatures") 
#Y0 = array(data[1:2,sample(1:dim(data)[2],1),],c(2,1,1))
#Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = 1)

## Not run
# Fit Non Homogeneous MS-AR models - univariate time series
#data(lynx)
#T = length(lynx)
#data = array(log10(lynx),c(T,1,1))
#theta.init = init.theta.MSAR(data,M=2,order=2,label="HH")
#mod.lynx.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=200)
#regimes.plot.MSAR(mod.lynx.hh,data,ylab="Captures number")
## End (not run)
}
