\name{simule.nh.MSAR.VM}
\alias{simule.nh.MSAR.VM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of (non) homogeneous Markov Stiwtching autoregressive models von Mises innovations
}
\description{
simule.nh.MSAR.VM simulates realisations of (non) homogeneous Markov Switching autoregressive models with  von Mises innovations}
\usage{
simule.nh.MSAR.VM(theta, Y0, T, N.samples = 1, covar.emis = NULL, covar.trans = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
list of class MSAR including model parameters and a description of the model. See init.theta.MSAR.VM for more details. 
}
  \item{Y0}{
  	Initial value. Array of dimension order*N.samples*d with order the AR order, N.samples the number of samples to be simulated and d the dimension of the considered data.
%%     ~~Describe \code{Y0} here~~
}
  \item{T}{Length of each realisation to be simulated
%%     ~~Describe \code{T} here~~
}
  \item{N.samples}{ number of samples to be simulated
%%     ~~Describe \code{N.samples} here~~
}
  \item{covar.emis}{ emission covariate or lag for non homogeneous models. Lag is used if the covariate is the lagged time series.
%%     ~~Describe \code{covar.emis} here~~
}
  \item{covar.trans}{transition covariate or lag for non homogeneous models. Lag is used if the covariate is the lagged time series.
%%     ~~Describe \code{covar.trans} here~~
}
}

\value{List including
	\item{..$Y}{simulated observation time series}
	\item{..$S}{simulated Markov chain}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ailliot P., Bessac J., Monbet V., P\`ene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI. 
}
\author{
Val\'erie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR.VM, init.theta.MSAR.VM}
\examples{
	##Not run
#data(WindDir)
#T = dim(WindDir)[1]
#N.samples = dim(WindDir)[2]
#Y = array(WindDir,c(T,N.samples,1))
# von Mises homogeneous MSAR
#M = 2
#order = 1
#theta.init = init.theta.MSAR.VM(Y,M=M,order=order,label="HH")
#polar.hh = fit.MSAR.VM(Y,theta.init,MaxIter=50,verbose=TRUE,eps=1e-8)

#K.sim = 1
#Y0 = array(Y[1:2,sample(1:N.samples,K.sim,replace=T),],c(2,K.sim,1))
#sim.dir = simule.nh.MSAR.VM(polar.hh$theta,Y0=Y0,T,N.samples=K.sim)

## Not run
#theta.init$mu = polar.hh$theta$mu
# theta.init$kappa = polar.hh$theta$kappa+1i*0 # kappa complex
# theta.init$prior = polar.hh$theta$prior
# theta.init$transmat = polar.hh$theta$transmat
# polar.hh.c = fit.MSAR.VM(Y,theta.init,MaxIter=50,verbose=TRUE,eps=1e-8)
	
# theta.init = init.theta.MSAR.VM(Y,M=M,order=order,label="NH",ncov=1,nh.transitions="VM")
# theta.init$mu = polar.hh.c$theta$mu
# theta.init$kappa = polar.hh.c$theta$kappa # kappa complex
# theta.init$prior = polar.hh.c$theta$prior
# theta.init$transmat = polar.hh.c$theta$transmat
# theta.init$par.trans =  matrix(c(polar.hh.c$theta$mu,.1*matrix(1,M,1)),M,2)+1i
#Y.tmp = array(Y[2:T,,],c(T-1,N.samples,1))
#Z = array(Y[1:(T-1),,],c(T-1,N.samples,1))
# polar.nh.c = fit.MSAR.VM(Y.tmp,theta.init,MaxIter=1,verbose=T,eps=1e-8,covar.trans=Z)
#K.sim = 100
#Y0 = array(Y[1:2,sample(1:N.samples,K.sim,replace=T),],c(2,K.sim,1))
#sim.dir = simule.nh.MSAR.VM(polar.nh.c$theta,Y0=Y0,T,N.samples=K.sim,covar.trans=1)


}
\keyword{ Simulation }
\keyword{ MSAR }
\keyword{ von Mises }% __ONLY ONE__ keyword per line
