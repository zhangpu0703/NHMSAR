\name{log_dens_Von_Mises}
\alias{log_dens_Von_Mises}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{von Mises log likelihood.}
\description{von Mises log likelihood.}
\usage{log_dens_Von_Mises(x, m, k)}
\arguments{
  \item{x}{vector of data}
  \item{m}{location parameter}
  \item{k}{dispersion parameter}
}
\details{
Log-likelihood of von Mises distribution with density 
\deqn{\frac{exp(k cos(x-m))}{2 \pi I_0(k)}}
where \eqn{I_0} is the modified Bessel function of order 0. }
\value{log likelihood}
\references{
Mardia, K.; Jupp, P. E. (1999). Directional Statistics. Wiley.}
\author{Valerie Monbet, valerie.monbet at univ-rennes1.fr}

\seealso{circular package}


\keyword{ von Mises }
\keyword{log-likelihood }