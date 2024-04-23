\name{lad}
\alias{lad}
\title{Least absolute deviations regression}
\description{
  This function is used to fit linear models considering Laplace errors.
}
\usage{
lad(formula, data, subset, na.action, method = "BR", tol = 1e-7, maxiter = 200,
  x = FALSE, y = FALSE, contrasts = NULL)
}
\arguments{
  \item{formula}{ an object of class \code{"formula"}: a symbolic description of
    the model to be fitted.}
  \item{data}{ an optional data frame containing the variables in the model. If
    not found in \code{data}, the variables are taken from \code{environment(formula)},
    typically the environment from which \code{lad} is called.}
  \item{subset}{ an optional expression indicating the subset of the rows of
    data that should be used in the fit.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{method}{ character string specifying the fitting method to be used; the options
    are \code{"BR"} Barrodale and Roberts' method (the default) and \code{"EM"} for an EM
    algorithm using IRLS.}
  \item{tol}{ the relative tolerance for the iterative algorithm. Default is \code{tol = 1e-7}.}
  \item{maxiter}{ The maximum number of iterations for the \code{EM} method. Default to 200.}
  \item{x, y}{ logicals. If \code{TRUE} the corresponding components of the fit (the model matrix, 
    the response) are returned.}
  \item{contrasts}{ an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
}
\value{
  An object of class \code{lad} representing the linear model fit. Generic
  function \code{print}, show the results of the fit.

  The functions \code{print} and \code{summary} are used to obtain and print a summary
  of the results. The generic accessor functions \code{coefficients}, \code{fitted.values}
  and \code{residuals} extract various useful features of the value returned by \code{lad}.
}
\author{The design was inspired by the R function \code{\link{lm}}.}
\references{
  Barrodale, I., Roberts, F.D.K. (1974).
  Solution of an overdetermined system of equations in the L1 norm.
  \emph{Communications of the ACM} \bold{17}, 319-320.

  Phillips, R.F. (2002).
  Least absolute deviations estimation via the EM algorithm.
  \emph{Statistics and Computing} \bold{12}, 281-285.
}
\examples{
fm <- lad(stack.loss ~ ., data = stackloss, method = "BR")
summary(fm)
}
\keyword{regression}
