\name{Laplace}
\alias{Laplace}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
\title{The Laplace distribution}
\description{
  Density, distribution function, quantile function and random generation for the
  Laplace distribution with location parameter \code{location} and scale parameter
  \code{scale}.
}
\usage{
dlaplace(x, location = 0, scale = 1, log = FALSE)
plaplace(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qlaplace(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rlaplace(n, location = 0, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{location, scale}{location and scale parameters. Scale must be positive.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]},
    otherwise, \eqn{P[X > x]}.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is taken
    to be the number required.}
}
\value{
  \code{dlaplace}, \code{plaplace}, and \code{qlaplace} are respectively the density,
  distribution function and quantile function of the Laplace distribution. \code{rlaplace}
  generates random deviates from the Laplace.

  The length of the result is determined by \code{n} for \code{rlaplace}, and is
  the maximum of the lengths of the numerical parameters for the other functions.
}
\details{
  If \code{location} or \code{scale} are not specified, they assume the default
  values of \code{0} and \code{1} respectively.

  The Laplace distribution with location \eqn{\mu} and scale \eqn{\phi} has density
  \deqn{
    f(x) =
    \frac{1}{\sqrt{2}\phi} \exp(-\sqrt{2}|x-\mu|/\phi)}{
    f(x) = 1/(\phi\sqrt2) exp(-\sqrt2|x - \mu|/\phi)
  }
}
\author{Felipe Osorio and Tymoteusz Wolodzko}
\references{
  Kotz, S., Kozubowski, T.J., and Podgorski, K. (2001).
  \emph{The Laplace Distributions and Generalizations}.
  Birkhauser, Boston.

  Phillips, R.F. (2002).
  Least absolute deviations estimation via the EM algorithm.
  \emph{Statistics and Computing} \bold{12}, 281-285.
}
\seealso{
  \link{Distributions} for other standard distributions and \code{\link{rmLaplace}}
  for the random generation from the multivariate Laplace distribution.
}
\examples{
x <- rlaplace(1000)
## Q-Q plot for Laplace data against true theoretical distribution:
qqplot(qlaplace(ppoints(1000)), x, main = "Laplace Q-Q plot",
  xlab = "Theoretical quantiles", ylab = "Sample quantiles")
abline(c(0,1), col = "red", lwd = 2)
}
\keyword{distribution}
