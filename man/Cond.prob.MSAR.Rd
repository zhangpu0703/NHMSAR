\name{Cond.prob.MSAR}
\alias{Cond.prob.MSAR}
\title{
Conditional probabilities for (non) homogeneous MSAR models}
\description{
Computes, for each time \eqn{t}, the conditional probabilities for MSAR models \eqn{P(Y_t|y_{1:(t-1)},y_{(t+1):T})} where  \eqn{Y} is the observed process and \eqn{y} the observed time series.}
\usage{
Cond.prob.MSAR(data, theta, yrange = NULL, covar.emis = NULL, covar.trans = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed time series, array of dimension T*N.samples*d 
}
  \item{theta}{
object of class MSAR including the model's parameter and description. See init.theta.MSAR for more details. }
  \item{yrange}{
values at which to compute the conditional probabilities
}
  \item{covar.emis}{
emission covariate if any.}
  \item{covar.trans}{
transition covariate if any.}
}

\value{a list including
	\item{..$yrange}{values at which the conditional probabilities are computed}
	\item{..$prob}{conditional probabilities for each time t and each values of yrange}
	\item{..$Yhat}{mode of the conditinal distribution for each time t}
}
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}
\seealso{
predict.MSAR}
\examples{
data(lynx)
data = array(log10(lynx),c(length(lynx),1,1))
T = length(data)
theta.init = init.theta.MSAR(data,M=2,order=2,label="HH")
mod.lynx.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=200)
ex = 100:114
lex = length(ex)
tps = (1821:1934)[ex] 
CP = Cond.prob.MSAR(array(data[ex,,],c(lex,1,1)), mod.lynx.hh$theta)
par(mfrow=c(2,1))
plot(tps,data[ex],typ="l",main="Homogeneous MSAR model",xlab="Time",ylab="Captured")
lines(tps,CP$Yhat,col="red")
alpha = .05
IC.emp = matrix(0,2,lex)
for (k in 1:lex) {
	tmp = cumsum(CP$prob[,k,])/sum(CP$prob[,k,])
	IC.emp[1,k] = CP$yrange[max(c(which(tmp<alpha/2),1))]
	IC.emp[2,k] = CP$yrange[min(max(which(tmp<(1-alpha/2))),length(CP$yrange))]
}
lines(tps,IC.emp[1,],lty=2,col="red")
lines(tps,IC.emp[2,],lty=2,col="red")

## Not run
#order = 2
#theta.init = init.theta.MSAR(data,M=2,order=2,label="NH",nh.transitions="logistic")
#theta.init$A0 = mod.lynx.hh$theta$A0
#theta.init$A = mod.lynx.hh$theta$A
#theta.init$sigma = mod.lynx.hh$theta$sigma
#theta.init$transmat = mod.lynx.hh$theta$transmat
#theta.init$prior = mod.lynx.hh$theta$prior
Y = array(data[2:T,,],c(T-1,1,1))
Z = array(data[1:(T-1),,],c(T-1,1,1))
#mod.lynx = fit.MSAR(array(Y,theta.init,covar.trans=Z)
Y.ex = array(data[ex,,],c(lex,1,1))
Z.ex = array(data[ex-1,,],c(lex,1,1))
#CPnh = Cond.prob.MSAR(Y.ex,mod.lynx$theta,covar.trans = Z.ex)
#
#plot(tps,data[ex],typ="l",main="Non Homogeneous MSAR model",xlab="Time",ylab="Captured")
#lines(tps,CPnh$Yhat,col="red")
#IC.emp = matrix(0,2,lex)
#for (k in 1:lex) {
#	tmp = cumsum(CPnh$prob[,k,])/sum(CPnh$prob[,k,])
#	IC.emp[1,k] = CPnh$yrange[max(c(which(tmp<alpha/2),1))]
#	IC.emp[2,k] = CPnh$yrange[min(max(which(tmp<(1-alpha/2))),length(CP$yrange))]
#}
#lines(tps,IC.emp[1,],lty=2,col="red")
#lines(tps,IC.emp[2,],lty=2,col="red")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MSAR }
\keyword{ Conditionnal probabilities }% __ONLY ONE__ keyword per line
