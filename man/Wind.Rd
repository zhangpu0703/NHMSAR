\name{Wind}
\alias{U}
\docType{data}
\title{
Winter wind data at 18 locations offshore of France}
\description{
Wind intensity at 18 locations offshore of France for months january and february. 32 years of data.  Time step is 6 hours. }
\usage{data(meteo.data)}
\format{
An array of dimension 248*32*18
  \describe{
    \item{\code{U}}{wind intensity}
   }
}

\source{
ERA-Interim}
\references{
Bessac, J., Ailliot, P., & Monbet, V. (2013). Gaussian linear state-space model for wind fields in the North-East Atlantic. arXiv preprint arXiv:1312.5530.}
\examples{
data(Wind)
}
\keyword{datasets}
