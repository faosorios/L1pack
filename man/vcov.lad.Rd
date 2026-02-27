\name{vcov.lad}
\alias{vcov.lad}
\title{Calculate variance-covariance matrix from \code{lad} models}
\description{
  Returns the variance-covariance matrix of the main parameters of a fitted model
  for \code{lad} objects. The \dQuote{main} parameters of model correspond to those
  returned by \code{\link{coef}}, and typically do not contain the nuisance scale
  parameter.
}
\usage{
\S3method{vcov}{lad}(object, level = 0.95, \dots)
}
\arguments{
  \item{object}{ an object representing a fitted model.}
  \item{level}{ the confidence level required for estimating the \emph{nuisance} 
    parameter involved in the asymptotic covariance matrix.}
  \item{\dots}{ additional arguments for method functions.}
}
\value{
  A matrix of the estimated covariances between the parameter estimates in the
  linear regression model. This should have row and column names corresponding
  to the parameter names given by the \code{\link{coef}} method. The nuisance 
  parameter required for the evaluation of the asympotic covariance matrix is 
  estimated using the method proposed by McKean and Schrader (1987).
}
\references{
  Gonin, R., Money, A.H. (1989).
  \emph{Nonlinear Lp-norm Estimation}.
  Marcel Dekker, New York.

  McKean, J.W., Schrader, R.M. (1987). Least absolute errors analysis of variance.
  In: Y. Dodge (Ed.), \emph{Statistical data analysis - Based on the L1-norm and related methods}. 
  North Holland, Amsterdam, pp. 297-305.
}
\keyword{models}
