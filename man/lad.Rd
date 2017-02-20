\name{lad}
\alias{lad}
\title{Least absolute deviations regression}
\description{
  This function is used to fit linear models considering Laplace errors.
}
\usage{
lad(formula, data, method = c("BR", "EM"), subset, na.action,
  control, model = TRUE, x = FALSE, y = FALSE, contrasts = NULL)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description of 
  the model to be fitted.}
  \item{data}{an optional data frame containing the variables in the model. If 
  not found in \code{data}, the variables are taken from \code{environment(formula)}, 
  typically the environment from which \code{lad} is called.}
  \item{method}{character string specifying the algorithm to use. The default 
  algorithm is the Barrodale and Roberts algorithm \code{method = "BR"}. Other 
  possible value is \code{method = "EM"} for an EM algorithm using IRLS.}
  \item{subset}{an optional expression indicating the subset of the rows of 
  data that should be used in the fit.}
  \item{na.action}{a function that indicates what should happen when the data contain NAs.}
  \item{control}{a list of control values for the estimation algorithm to replace 
  the default values returned by the function \code{\link{l1pack.control}}.}
  \item{model, x, y}{logicals. If \code{TRUE} the corresponding components of 
  the fit (the model frame, the model matrix, the response) are returned.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
}
\value{
  an object of class \code{lad} representing the linear model fit. Generic 
  function \code{print}, show the results of the fit.
}
\author{The design was inspired by the R function \code{\link{lm}}.}
\references{
  Barrodale, I., and Roberts, F.D.K. (1974). 
  Solution of an overdetermined system of equations in the L1 norm. 
  \emph{Communications of the ACM} \bold{17}, 319-320.

  Phillips, R.F. (2002).
  Least absolute deviations estimation via the EM algorithm.
  \emph{Statistics and Computing} \bold{12}, 281-285.
}
\examples{
lad(stack.loss ~ ., data = stackloss, method = "EM")
}
\keyword{regression}
