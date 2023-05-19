\name{lad.fit}
\alias{lad.fit}
\title{Fitter functions for least absolute deviation (LAD) regression}
\usage{
lad.fit(x, y, method = "BR", tol = 1e-7, maxiter = 200)
}
\description{
  This function is a \emph{switcher} among various numerical fitting functions
  (\code{\link{lad.fit.BR}}, and \code{\link{lad.fit.EM}}). The argument \code{method}
  does the switching: \code{"BR"} for \code{\link{lad.fit.BR}}, etc. This should
  usually \emph{not} be used directly unless by experienced users.
}
\arguments{
  \item{x}{ design matrix of dimension \eqn{n\times p}.}
  \item{y}{ vector of observations of length \eqn{n}.}
  \item{method}{ currently, methods \code{"BR"} (default), and \code{"EM"} are supported.}
  \item{tol}{ the relative tolerance for the iterative algorithm. Default is \code{tol = 1e-7}.}
  \item{maxiter}{ The maximum number of iterations for the \code{EM} method. Default to 200.}
}
\value{
  a \code{\link{list}} with components:
  \item{coefficients}{ a named vector of coefficients.}
  \item{scale}{ final scale estimate of the random error.}
  \item{residuals}{ the residuals, that is response minus fitted values.}
  \item{fitted.values}{ the fitted values.}
  \item{SAD}{ the sum of absolute deviations.}
  \item{weights}{ estimated \code{EM} weights.}
  \item{basic}{ basic observations, that is observations with zero residuals.}
  \item{logLik}{ the log-likelihood at convergence.}
}
\seealso{
  \code{\link{lad.fit.BR}}, \code{\link{lad.fit.EM}}.
}
\examples{
x <- cbind(1, stack.x)
fm <- lad.fit(x, stack.loss, method = "BR")
fm
}
\keyword{regression}
