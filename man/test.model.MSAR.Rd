\name{test.model.MSAR}
\alias{test.model.MSAR}
\title{
Performs bootstrap statistical tests to validate MSAR models. }
\description{
Performs bootstrap statistical tests to validate MSAR models. Marginal distribution, auto correlation function and up-crossings are considered. For each of them the tests statistic computed from observations is compared to the distribution of the satistics corresponding to the MSAR model.}
\usage{
test.model.MSAR(data,simu,lag=NULL,id=1,u=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*d }
  \item{simu}{
  simulated time series, array of dimension T*N.sim*d.
N.sim have to be  K*N.samples with K large enough (for instance, K=100) }
 \item{lag}{maximum lag for auto-correlation functions. }
 \item{id}{considered component. It is usefull when data is multivariate. }
 \item{u}{considered levels for up crossings}
 }

\value{Returns a list including
	\item{StaDist}{statistics of marginal distributions, based on Smirnov like statistics}
	\item{..$dd}{test statistic} 
	\item{..$q.dd}{quantiles .05 and .95 of the distribution of the test statistic under the null  hypothesis } 
	\item{..$p.value}{p value} 
	\item{Cor}{statistics of correlation functions}
	\item{..$dd}{test statistic} 
	\item{..$q.dd}{quantiles .05 and .95 of the distribution of the test statistic under the null  hypothesis } 
	\item{..$p.value}{p value} 
	\item{ENu}{statistics of  intensity of up crossings }
	\item{..$dd}{test statistic} 
	\item{..$q.dd}{quantiles .05 and .95 of the distribution of the test statistic under the null  hypothesis } 
	\item{..$p.value}{p value} 
	\item{AD}{statistics of marginal distributions, based on Anderson Darling statistics}
	\item{..$dd}{test statistic} 
	\item{..$q.dd}{quantiles .05 and .95 of the distribution of the test statistic under the null  hypothesis } 
	\item{..$p.value}{p value} 

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\details{ Test statistics
	Marginal distribution: {\deqn{ S = \int_{-\infty}^{\infty} \left| F_n(x)-F(x) \right| dx} } 
	
	Marginal distribution, based on Anderson Darling statistic: {\deqn{ S = \int_{-\infty}^{\infty} \left| \frac{F_n(x)-F(x)}{F(x)(1-F(x))} \right| dx} } 
	
		Correlation function: {\deqn{ S = \int_0^L\left|C_n(l)-C(l)\right|dl} } 

	Number of up crossings: {\deqn{ S = \int_{-\infty}^{\infty}\left|E_n(N_u)-E(N_u)\right|du}}
}
%\references{
%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}




%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MSAR model validation }
\seealso{valid_all, test.model.MSAR}