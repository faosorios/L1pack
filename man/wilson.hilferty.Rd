\name{WH.Laplace}
\alias{WH.Laplace}
\title{Wilson-Hilferty transformation}
\description{
  Returns the Wilson-Hilferty transformation for multivariate Laplace deviates.
}
\usage{
WH.Laplace(x, center, Scatter)
}
\arguments{
  \item{x}{ object of class \code{'LaplaceFit'} from which is extracted the estimated Mahalanobis distances of the fitted model.
  Also \code{x} can be a vector or matrix of data with, say, \eqn{p} columns.}
  \item{center}{ mean vector of the distribution or data vector of length \eqn{p}. Not required if \code{x} have class \code{'LaplaceFit'}.}
  \item{Scatter}{ Scatter matrix (\eqn{p} by \eqn{p}) of the distribution. Not required if \code{x} have class \code{'LaplaceFit'}.}
}
\details{
  Let \eqn{T = D/(2p)}{T = D/(2p)} be a Gamma distributed random variable, where \eqn{D^2} 
  denotes the squared Mahalanobis distance defined as
  \deqn{D^2 = (\bold{x} - \bold{\mu})^T \bold{\Sigma}^{-1} (\bold{x} - \bold{\mu}).}

  Thus, the Wilson-Hilferty transformation is given by
  \deqn{z = \frac{T^{1/3} - (1 - \frac{1}{9p})}{(\frac{1}{9p})^{1/2}}%
  }{z = (T^1/3 - (1 - 1/(9p))) / (1/(9p))^1/2}
  and \eqn{z} is approximately distributed as a standard normal distribution. This is useful, 
  for instance, in the construction of QQ-plots.
}
\references{
  Osorio, F., Galea, M., Henriquez, C., Arellano-Valle, R. (2023).
  Addressing non-normality in multivariate analysis using the \emph{t}-distribution. 
  \emph{AStA Advances in Statistical Analysis} \bold{107}, 785-813.

  Terrell, G.R. (2003).
  The Wilson-Hilferty transformation is locally saddlepoint.
  \emph{Biometrika} \bold{90}, 445-453.

  Wilson, E.B., Hilferty, M.M. (1931).
  The distribution of chi-square.
  \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{17}, 684-688.
}
\examples{
Scatter <- matrix(c(1,.5,.5,1), ncol = 2)
Scatter

# generate the sample
y <- rmLaplace(n = 500, Scatter = Scatter)
fit <- LaplaceFit(y)
z <- WH.Laplace(fit)
par(pty = "s")
qqnorm(z, main = "Transformed distances QQ-plot")
abline(c(0,1), col = "red", lwd = 2)
}
\keyword{multivariate}
