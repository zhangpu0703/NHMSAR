\name{simule_MC}
\alias{simule_MC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates Markov chain of length T}
\description{
Simulates Markov chain of length T, given a transition matrix and a prior distribution. }
\usage{
simule_MC(transmat, prior, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transmat}{
transition matrix}
  \item{prior}{
prior distribution}
  \item{T}{
simulation length}
}

\value{
\item{X}{Markov chain sequence} }
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
simule_MC.nh, simule.nh.MSAR}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
\keyword{ Markov chain }% __ONLY ONE__ keyword per line
