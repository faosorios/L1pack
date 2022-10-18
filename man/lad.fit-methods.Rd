\name{lad.fit-methods}
\alias{lad.fit.BR}
\alias{lad.fit.EM}
\title{Fit a least absolute deviation (LAD) regression model}
\usage{
lad.fit.BR(x, y, tol = 1e-7)
lad.fit.EM(x, y, tol = 1e-7, maxiter = 500)
}
\description{
  Fits a linear model using LAD methods, returning the bare minimum computations.
}
\arguments{
  \item{x, y}{ numeric vectors or matrices for the predictors and the response in
  a linear model. Typically, but not necessarily, \code{x} will be constructed by
  one of the fitting functions.}
  \item{tol}{ the relative tolerance for the iterative algorithm. Default is \code{tol = 1e-7}.}
  \item{maxiter}{ The maximum number of iterations for the \code{EM} method. Default to 500.}
}
\value{
  The bare bones of a \code{lad} object: the coefficients, residuals, fitted values,
  and some information used by \code{summary.lad}.
}
\seealso{
  \code{\link{lad}}, \code{\link{lad.fit}}, \code{\link{lm}}
}
\examples{
x <- cbind(1, stack.x)
z <- lad.fit.BR(x, stack.loss)
z
}
\keyword{regression}
