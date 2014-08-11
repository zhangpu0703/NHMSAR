\name{WindDir}
\alias{WindDir}
\docType{data}
\title{
January wind direction at Ouessant}
\description{
Wind direction at Ouessant. 49 independant january month (one per column). Time step is 6 hours. }
\usage{data(meteo.data)}
\format{
A matrix of dimension 124*32
  \describe{
    \item{\code{WindDir}}{wind direction}
   }
}

\source{
ERA-Interim}
\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.}
\examples{
data(WindDir)
}
\keyword{datasets}
