\name{mLaplace}
\alias{mLaplace}
\alias{dmLaplace}
\alias{rmLaplace}

\title{Multivariate Laplace distribution}
\description{
  These functions provide the density and random number generation from the multivariate 
  Laplace distribution.
}
\usage{
dmLaplace(x, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)), log = FALSE)
rmLaplace(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
}

\arguments{
  \item{x}{ vector or matrix of data. }
  \item{n}{ the number of samples requested.}
  \item{center}{ a \eqn{k\times 1}{k x 1} vector giving the locations.}
  \item{Scatter}{ a \eqn{k \times k}{k x k} positive-definite dispersion matrix.}
  \item{log}{ logical; if TRUE, the logarithm of the density function is returned.}
}

\value{
  If \code{x} is a matrix with \eqn{n} rows, then \code{dmLaplace} returns a 
  \eqn{n\times 1}{n x 1} vector considering each row of \code{x} as a copy from 
  the multivariate Laplace.

  If \code{n = 1}, then \code{rmLaplace} returns a vector of the same length as 
  \code{center}, otherwise a matrix of \code{n} rows of random vectors.
}

\details{
  The multivariate Laplace distribution is a multidimensional extension of the univariate 
  symmetric Laplace distribution. There are multiple forms of the multivariate Laplace 
  distribution. Here, a particular case of the multivarite power exponential distribution 
  introduced by Gomez et al. (1998) is considered.

  The multivariate Laplace distribution with location \eqn{\bold{\mu}} = \code{center} 
  and \eqn{\bold{\Sigma}} = \code{Scatter} has density
  \deqn{
    f(\bold{x}) =
    \frac{\Gamma(k/2)}{\pi^{k/2}\Gamma(k)2^{k+1}}|\bold{\Sigma}|^{-1/2}
    \exp\left\{-\frac{1}{2}[(\bold{x} - \bold{\mu})^T \bold{\Sigma}^{-1} (\bold{x} 
    - \bold{\mu})]^{1/2}\right\}.
  }
  The function \code{rmLaplace} is an interface to C routines, which make calls to
  subroutines from LAPACK. The matrix decomposition is internally done using
  the Cholesky decomposition. If \code{Scatter} is not non-negative definite then
  there will be a warning message.
}

\references{
  Fang, K.T., Kotz, S., Ng, K.W. (1990).
  \emph{Symmetric Multivariate and Related Distributions}.
  Chapman & Hall, London.

  Gomez, E., Gomez-Villegas, M.A., Marin, J.M. (1998).
  A multivariate generalization of the power exponential family of distributions.
  \emph{Communications in Statistics - Theory and Methods} \bold{27}, 589-600.
  \doi{10.1080/03610929808832115}
}

\examples{
# dispersion parameters
Scatter <- matrix(c(1,.5,.5,1), ncol = 2)
Scatter

# generate the sample
y <- rmLaplace(n = 2000, Scatter = Scatter)

# scatterplot of a random bivariate Laplace sample with center
# vector zero and scale matrix 'Scatter'
par(pty = "s")
plot(y, xlab = "", ylab = "")
title("bivariate Laplace sample", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
