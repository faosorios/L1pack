\name{l1pack.control}
\alias{l1pack.control}

\title{Set control parameters}
\description{
  Allows users to set parameters for \code{\link{lad}}.
}
\usage{
l1pack.control(maxIter = 2000, tolerance = 1e-9)
}

\arguments{
  \item{maxIter}{maximum number of iterations. The default is 2000.}
  \item{tolerance}{the relative tolerance in the iterative algorithm.}
}

\value{
  A list of control arguments to be used in a call to \code{\link{lad}}.

  A call to \code{l1pack.control} can be used directly in the control argument 
  of a call to \code{\link{lad}}. 
}
\author{Felipe Osorio.}
\examples{
ctrl <- l1pack.control(maxIter = 50, tol = 1e-07)
lad(stack.loss ~ ., data = stackloss, control = ctrl)
}
\keyword{regression}
\keyword{models}
