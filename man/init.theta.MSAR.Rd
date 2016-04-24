\name{init.theta.MSAR (NH-MSAR)}
\alias{init.theta.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialisation function for MSAR model fitting}
\description{
Initialization before fitting (non) homogeneous Markov switching autoregressive models by EM algorithm. Non homogeneity may be introduce at the intercept level or in the probability transitions. The link functions are defined here.}
\usage{
init.theta.MSAR(data, ..., M, order, regime_names = NULL, nh.emissions = NULL, 
nh.transitions = NULL, label = NULL, ncov.emis = 0, ncov.trans = 0,cl.init="mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d with 
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension}
 
  \item{M}{
number of regimes}
  \item{order}{
order of AR processes}
  \item{label}{
"HH" (default) for homogeneous MS AR model \\
"HN" for non homogeneous emissions \\
"NH" for non homogeneous transitions  \\
"NN" for non homogeneous emissions and non homogeneous transitions }
  \item{regime_names}{
(optional) regime's names may be chosen}
  \item{nh.emissions}{
link function for non homogeneous emissions. If nh.emissions="linear" (default) linear link is used. If you define an other function it should follow the sample nh.emissions <- function(covar,par.emis) with par.emis of dimension M by ncov.emis+1.   }
  \item{nh.transitions}{
link function for non homogeneous transitions. If nh.transitions="gauss" (default) gaussian link is used. If M=2, "logistic" may be chosen. If you define an other function it should follow the sample nh.transitions <- function(covar,par.trans,transma) with par.emis of dimension M by ncov.trans+1.  }

  \item{ncov.emis}{
number of covariates in HN model}
  \item{ncov.trans}{
number of covariates in NH model
}
\item{cl.init}{allows to choose the initialization method. }
 \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The default implemented link function for non homogneneous intercept is the linear function \deqn{A0_t^{(x)} = \theta_{A0}^{(x)} Z(t)}
\eqn{\theta_{A0}^{(x)}} denotes a line vector here.
Other link functions can be defined using nh.emissions (see above).

The default implemented link function for non homogneneous transitions is the Gauss function. Transition from  \eqn{i} to \eqn{j} is defined as follows. 

\deqn{f(Z,\theta_Q,Q;i,j) = Q_{ij} \exp\left(-\frac{1}{2}\frac{(Z-\theta_Q^{(j)}(1))^2}{\theta_Q^{(j)}(2)}\right)} then \eqn{f} is normalized in order to define a stochastic matrix. 

When, only two regimes are considered, the logistic link can be used. Probability of staying in state \eqn{i} is defined as follows 
\deqn{f(Z,\theta_Q,Q;i,i) = \epsilon+(-2-\epsilon)/(1+\exp(\theta_Q^{(i)}(1)+\theta_Q^{(i)}[2:(d_Z+1)] Z)) }
\deqn{f(Z,\theta_Q,Q;i,j) = 1-f(Z,\theta_Q,Q;i,i)}
with \eqn{Z} the covariate and eqn{d_Z} its dimension (number of covariates)
}
\value{return a list of class MSAR including 
	
\item{theta}{parameter}
\item{..$transmat}{transition matrix}
\item{..$prior}{prior probabilities}
\item{..$A}{list including the autoregressive coefficients (or matrices)}
\item{..$A0}{intercepts}
\item{..$sigma}{variances of innovations}
\item{..$par.emis}{parameters of non homogeneous emissions} 
\item{..$par.trans}{parameters of non homogeneous transitions} 
\item{label}{model's label}
}
\references{
Ailliot, Monbet}
\author{
Val\'erie Monbet, valerie.monbet at univ-rennes1.fr}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{fit.MSAR}
\examples{

data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]

# Fit Homogeneous MS-AR models
M = 2
order = 2
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH") 
mod.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=10)
regimes.plot.MSAR(mod.hh,data,ylab="temperatures") 

## Not run
# Fit Non Homogeneous MS-AR models
#theta.init = init.theta.MSAR(data,M=M,order=order,label="NH",nh.transitions="gauss")
#attributes(theta.init)
#mod.nh = fit.MSAR(array(data[2:T,,],c(T-1,N.samples,1)),theta.init,verbose=TRUE,MaxIter=50,
#covar.trans=array(data[1:(T-1),,],c(T-1,N.samples,1)))
#regimes.plot.MSAR(mod.nh,data,ex=40,ylab="temperature (deg. C)")

## Not run
# Fit Non Homogeneous MS-AR models to lynx data
#data(lynx)
#data = array(lynx,c(length(lynx),1,1))
#theta.init = init.theta.MSAR(data,M=2,order=2,label="NH",nh.transitions="logistic")
#attributes(theta.init)
#mod.lynx = fit.MSAR(array(data[2:T,,],c(T-1,1,1)),theta.init,verbose=TRUE,MaxIter=200,
#covar.trans=array(data[1:(T-1),,],c(T-1,1,1)))
#regimes.plot.MSAR(mod.lynx,data,ylab="Captures number")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ initialisation }
\keyword{ MSAR model fitting }% __ONLY ONE__ keyword per line
