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
\S3method{vcov}{lad}(object, \dots)
}
\arguments{
  \item{object}{ an object representing a fitted model.}
  \item{\dots}{ additional arguments for method functions.}
}
\value{
  A matrix of the estimated covariances between the parameter estimates in the
  linear regression model. This should have row and column names corresponding
  to the parameter names given by the \code{\link{coef}} method.
}
\keyword{models}
