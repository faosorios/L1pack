\name{l1fit}
\alias{l1fit}
\title{Minimum absolute residual (L1) regression}
\description{
  Performs an L1 regression on a matrix of explanatory variables and a vector of responses.
}

\usage{
l1fit(x, y, intercept = TRUE, tolerance = 1e-07, print.it = TRUE)
}
\arguments{
  \item{x}{ vector or matrix of explanatory variables. Each row corresponds to an
    observation and each column to a variable. The number of rows of \code{x} should
    equal the number of data values in \code{y}, and there should be fewer columns
    than rows. Missing values are not allowed.}
  \item{y}{ numeric vector containing the response. Missing values are not allowed.}
  \item{intercept}{ logical flag. If \code{TRUE}, an intercept term is included in
    the regression model.}
  \item{tolerance}{ numerical value used to test for singularity in the regression.}
  \item{print.it}{ logical flag. If \code{TRUE}, then warnings about non-unique solutions
    and rank deficiency are given.}
}
\value{
  list defining the regression (compare with function \code{\link{lsfit}}).
  \item{coefficients}{ vector of coefficients.}
  \item{residuals}{ residuals from the fit.}
  \item{message}{ character strings stating whether a non-unique solution is possible, 
  or if the \code{x} matrix was found to be rank deficient.}
}
\details{%
  The Barrodale-Roberts algorithm, which is a specialized linear programming
  algorithm, is used.
}

\references{
  Barrodale, I., Roberts, F.D.K. (1973).
  An improved algorithm for discrete L1 linear approximations.
  \emph{SIAM Journal of Numerical Analysis} \bold{10}, 839-848.
  \doi{10.1137/0710069}

  Barrodale, I., Roberts, F.D.K. (1974).
  Solution of an overdetermined system of equations in the L1 norm.
  \emph{Communications of the ACM} \bold{17}, 319-320.
  \doi{10.1145/355616.361024}

  Bloomfield, P., Steiger, W.L. (1983).
  \emph{Least Absolute Deviations: Theory, Applications, and Algorithms.}
  Birkhauser, Boston, Mass.
}

\examples{
l1fit(stack.x, stack.loss)
}
\keyword{regression}
