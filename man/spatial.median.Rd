\name{spatial.median}
\alias{spatial.median}

\title{Computation of the generalized spatial median}
\description{
  Computation of the generalized spatial median estimator as defined by Rao (1988).
}
\usage{
spatial.median(x, data, subset, na.action, tol = 1e-6, maxiter = 200)
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
  A list with class \code{'spatial.median'} containing the following components:
  \item{call}{  a list containing an image of the \code{spatial.median} call that produced the object.}
  \item{median}{  final estimate of the location vector.}
  \item{Scatter}{  final estimate of the scale matrix.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{innerIter}{  the total number of iterations used in the \bold{inner} iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the Kotz distribution.}
  \item{distances}{  estimated squared Mahalanobis distances.}
  Generic function \code{print} show the results of the fit.
}

\details{
  An interesting fact is that the generalized spatial median estimator proposed by Rao (1988) is 
  the maximum likelihood estimator under the Kotz-type distribution discussed by Naik and Plungpongpun (2006).
  The generalized spatial median estimators are defined as \eqn{\hat{\bold{\mu}}} and \eqn{\hat{\bold{\Sigma}}}
  which minimize
  \deqn{
    \frac{n}{2}\log|\bold{\Sigma}| + \sum\limits_{i=1}^n \sqrt{(\bold{x} - \bold{\mu})^T 
    \bold{\Sigma}^{-1} (\bold{x} - \bold{\mu})},
  }
  simultaneously with respect to \eqn{\bold{\mu}} and \eqn{\bold{\Sigma}}.
  
  The function \code{spatial.median} follows the iterative reweighting algorithm of Naik and Plungpongpun (2006).
}

\references{
  Naik, D.N., Plungpongpun, K. (2006). 
  A Kotz-type distribution for multivariate statistical inference. 
  In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (Eds) 
  \emph{Advances in Distribution Theory, Order Statistics, and Inference}. 
  Birkhauser Boston, pp. 111-124.

  Rao, C.R. (1988). 
  Methodology based on the L1-norm in statistical inference.
  \emph{Sankhya, Series A} \bold{50}, 289-313. 
}

\seealso{
\code{\link{cov}}, \code{\link{LaplaceFit}}
}

\examples{
z <- spatial.median(stack.x)
z
}
\keyword{models}
\keyword{multivariate}
