\name{Laplace}
\alias{Laplace}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
\title{The symmetric Laplace distribution}
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
  \item{x, q}{ vector of quantiles.}
  \item{p}{ vector of probabilities.}
  \item{n}{ number of observations.}
  \item{location}{ location parameter \eqn{\mu}, which is the mean.}
  \item{scale}{ scale parameter \eqn{\phi}. Scale must be positive.}
  \item{log, log.p}{ logical; if TRUE, probabilities \eqn{p} are given as \eqn{\log(p)}.}
  \item{lower.tail}{ logical; if TRUE (default), probabilities are \eqn{P(X \le x)},
    otherwise, \eqn{P(X > x)}.}
}
\value{
  \code{dlaplace}, \code{plaplace}, and \code{qlaplace} are respectively the density,
  distribution function and quantile function of the Laplace distribution. \code{rlaplace}
  generates random deviates drawn from the Laplace distribution, the length of the result 
  is determined by \code{n}.
}
\details{
  If \code{location} or \code{scale} are not specified, they assume the default
  values of \code{0} and \code{1} respectively.

  The Laplace distribution with location \eqn{\mu} and scale \eqn{\phi} has density
  \deqn{
    f(x) =
    \frac{1}{\sqrt{2}\phi} \exp(-\sqrt{2}|x-\mu|/\phi),}{
    f(x) = 1/(\phi\sqrt2) exp(-\sqrt2|x - \mu|/\phi),
  }
  where \eqn{-\infty < y < \infty}, \eqn{-\infty < \mu < \infty} and \eqn{\phi > 0}. 
  The mean is \eqn{\mu} and the variance is \eqn{\phi^2}.

  The cumulative distribution function, assumes the form
  \deqn{ 
    F(x) = 
    \left\{\begin{array}{ll}
    \frac{1}{2} \exp(\sqrt{2}(x - \mu)/\phi) & x < \mu, \\
    1 - \frac{1}{2} \exp(-\sqrt{2}(x - \mu)/\phi) & x \geq \mu.
    \end{array}\right.}{
    F(x) = [if x < mu:] 1/2 * exp(sqrt(2)(x-\mu)/\phi)
           [else:] 1 - 1/2 * exp(-sqrt(2)(x-\mu)/\phi)
  }

  The quantile function, is given by 
  \deqn{ 
    F^{-1}(p) = \left\{\begin{array}{ll}
    \mu + \frac{\phi}{\sqrt{2}} \log(2p) & p < 0.5, \\
    \mu - \frac{\phi}{\sqrt{2}} \log(2(1-p)) & p \geq 0.5.
    \end{array}\right.}{
    F^-1(p) = [if p < 0.5:] \mu + (\phi/sqrt(2)) * log(2*p)
              [else:] \mu - (\phi/sqrt(2)) * log(2*(1-p))
  }
}
\author{Felipe Osorio and Tymoteusz Wolodzko}
\references{
  Kotz, S., Kozubowski, T.J., Podgorski, K. (2001).
  \emph{The Laplace Distributions and Generalizations}.
  Birkhauser, Boston.

  Krishnamoorthy, K. (2006). 
  \emph{Handbook of Statistical Distributions with Applications}, 2nd Ed. 
  Chapman & Hall, Boca Raton.
}
\seealso{
  \link{Distributions} for other standard distributions and \code{\link{rmLaplace}}
  for the random generation from the multivariate Laplace distribution.
}
\examples{
x <- rlaplace(1000)
## QQ-plot for Laplace data against true theoretical distribution:
qqplot(qlaplace(ppoints(1000)), x, main = "Laplace QQ-plot",
  xlab = "Theoretical quantiles", ylab = "Sample quantiles")
abline(c(0,1), col = "red", lwd = 2)
}
\keyword{distribution}
