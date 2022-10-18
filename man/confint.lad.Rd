\name{confint.lad}
\alias{confint.lad}
\title{Confidence intervals from \code{lad} models}
\description{
  Computes confidence intervals for one or more parameters in a fitted
  model corresponding to a \code{lad} object.
}
\usage{
\method{confint}{lad}(object, parm, level = 0.95, ...)
}
\arguments{
  \item{object}{ a fitted model object.}
  \item{parm}{ a specification of which parameters are to be given confidence intervals,
    either a vector of numbers or a vector of names. If missing, all parameters are
    considered.}
  \item{level}{ the confidence level required.}
  \item{\dots}{ additional argument(s) for methods.}
}
\value{
  A matrix (or vector) with columns giving lower and upper confidence limits for
  each parameter. These will be labelled as \code{(1-level)/2} and \code{1 - (1-level)/2}
  in \% (by default 2.5\% and 97.5\%).
}
\details{
  \code{confint} is a generic function. Confidence intervals associated to \code{lad}
  objects are asymptotic, and needs suitable \code{\link{coef}} and \code{\link{vcov}}
  methods to be available.
}
\seealso{
  \code{\link[MASS:confint]{confint.glm}} and
  \code{\link[MASS:confint]{confint.nls}} in package \CRANpkg{MASS}.
}
\examples{
fm <- lad(stack.loss ~ ., data = stackloss, method = "BR")
confint(fm) # based on asymptotic normality
}
\keyword{models}
