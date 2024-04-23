\name{LaplaceFit}
\alias{LaplaceFit}

\title{Estimation of mean and covariance using the multivariate Laplace distribution}
\description{
  Estimates the mean vector and covariance matrix assuming the data came from a multivariate
  Laplace distribution.
}
\usage{
LaplaceFit(x, data, subset, na.action, tol = 1e-6, maxiter = 200)
}
\arguments{
  \item{x}{  a formula or a numeric matrix or an object that can be coerced to a numeric matrix.}
  \item{data}{  an optional data frame (or similar: see \code{\link{model.frame}}), used only if
  \code{x} is a formula. By default the variables are taken from \code{environment(formula)}.}
  \item{subset}{ an optional expression indicating the subset of the rows of
    data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{tol}{ the relative tolerance in the iterative algorithm.}
  \item{maxiter}{ maximum number of iterations. The default is 200.}
}
\value{
  A list with class \code{'LaplaceFit'} containing the following components:
  \item{call}{  a list containing an image of the \code{LaplaceFit} call that produced the object.}
  \item{center}{  final estimate of the location vector.}
  \item{Scatter}{  final estimate of the scale matrix.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the Laplace distribution.}
  \item{distances}{  estimated squared Mahalanobis distances.}
  Generic function \code{print} show the results of the fit.
}
\references{
  Yavuz, F.G., Arslan, O. (2018).
  Linear mixed model with Laplace distribution (LLMM).
  \emph{Statistical Papers} \bold{59}, 271-289.
}

\seealso{
\code{\link{cov}}
}

\examples{
fit <- LaplaceFit(stack.x)
fit

# covariance matrix
p <- fit$dims[2]
Sigma <- (4 * (p + 1)) * fit$Scatter
Sigma
}
\keyword{models}
\keyword{multivariate}
