\name{rmLaplace}
\alias{rmLaplace}

\title{Multivariate Laplace Random Deviates}
\description{
  Random number generation from the multivariate Laplace distribution.
}
\usage{
rmLaplace(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{center}{ a vector giving the locations of each variable}
  \item{Scatter}{ a positive-definite dispersion matrix}
}

\value{
  If \code{n = 1} a vector of the same length as \code{center}, otherwise a
  matrix of \code{n} rows of random vectors.
}

\details{
  The function \code{rmLaplace} is an interface to C routines, which make calls to
  subroutines from LAPACK. The matrix decomposition is internally done using
  the Cholesky decomposition. If \code{Scatter} is not non-negative definite then
  there will be a warning message.
}

\references{
  Gomez, E., Gomez-Villegas, M.A., and Marin, J.M. (1998).
  A multivariate generalization of the power exponential family of distributions.
  \emph{Communications in Statistics - Theory and Methods} \bold{27}, 589-600.

  Kotz, S., Kozubowski, T.J., and Podgorski, K. (2001).
  \emph{The Laplace Distributions and Generalizations}.
  Birkhauser, Boston.
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
