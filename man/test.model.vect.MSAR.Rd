\name{test.model.vect.MSAR}
\alias{test.model.vect.MSAR}
\title{
Performs bootstrap statistical tests on covariance to validate MSVAR models. }
\description{
Performs bootstrap statistical on covariance to validate MSVAR models.}
\usage{
test.model.vect.MSAR(data,simu,lag=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*d }
  \item{simu}{
  simulated time series, array of dimension T*N.sim*d.
N.sim have to be  K*N.samples with K large enough (for instance, K=100) }
 \item{lag}{ to be considered (usefull for state space models)}
 }

\value{Returns a list including
	\item{Cvect}{statistics of covariance}
	\item{..$dd}{test statistic} 
	\item{..$q.dd}{quantiles .05 and .95 of the distribution of the test statistic underthe null  hypothesis } 
	\item{..$p.value}{p value} 

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\details{ Test statistics
	\deqn{ S = || C_n-C ||} }
%\references{
%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}




%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MSAR model validation }
\seealso{valid_all, test.model.MSAR}