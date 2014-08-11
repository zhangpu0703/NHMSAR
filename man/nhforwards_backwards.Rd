\name{nhforwards_backwards}
\alias{nhforwards_backwards}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward Backward for MSAR models with non homogeneous transitions }
\description{
 Computes the posterior (or smoothing) probabilities in an homogenenous HMM or MSAR model using the forwards backwards algo.
 'filter_only' is an optional argument (default: 0). If 1, we do filtering, if 0, smoothing.
}
\usage{
nhforwards_backwards(prior, transition, obslik, filter_only = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{
rior probabilities
PRIOR(I) = Pr(X(1) = I)}
  \item{transition}{
non homogeneous transitions, one transition matrix for each time}
  \item{obslik}{
emission probabilities
OBSLIK(I,t) = Pr(Y(t) | X(t)=I)}
  \item{filter_only}{
optional argument (default: 0). If TRUE, we do filtering, if FALSE, smoothing (default).
}
}

\value{ \item{..$gamma}{smoothing probabilities 
  	P(X(t)|Y(0),...,Y(T))}
  \item{..$xi }{two steps smoothing probabilities 
  	 	P(X(t),X(t+1)|Y(0),...,Y(T))}
  \item{..$loglik }{log likelihood}
  \item{..$M }{Number of regimes}
  \item{..$alpha}{intermediate component in the FB algorithm (forward)}
  \item{..$beta }{intermediate component in the FB algorithm (backward)}
}
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR, Estep.MSAR}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Forward Backward }
\keyword{ Em algorithm }% __ONLY ONE__ keyword per line
\keyword{ Smoothing probabilities }