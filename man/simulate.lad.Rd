\name{simulate.lad}
\alias{simulate.lad}
\title{Simulate Responses from \code{lad} Models}
\description{
  Simulate one or more responses from the distribution corresponding to a
  fitted \code{lad} object.
}
\usage{
\method{simulate}{lad}(object, nsim = 1, seed = NULL, \dots)
}
\arguments{
  \item{object}{an object representing a fitted model.}
  \item{nsim}{number of response vectors to simulate. Defaults to \code{1}.}
  \item{seed}{an object specifying if and how the random number generator should
    be initialized (\sQuote{seeded}). For the "lad" method, either \code{NULL}
    or an integer that will be used in a call to \code{set.seed} before simulating
    the response vectors. If set, the value is saved as the \code{"seed"} attribute
    of the returned value. The default, \code{NULL} will not change the random
    generator state, and return \code{\link{.Random.seed}} as the \code{"seed"}
    attribute, see \sQuote{Value}.
  }
  \item{\dots}{additional optional arguments.}
}
\value{
  For the \code{"lad"} method, the result is a data frame with an attribute
  \code{"seed"}. If argument \code{seed} is \code{NULL}, the attribute is the
  value of \code{\link{.Random.seed}} before the simulation was started.
}
\author{Tymoteusz Wolodzko and Felipe Osorio}
\examples{
fm <- lad(stack.loss ~ ., data = stackloss)
sm <- simulate(fm, nsim = 4)
}
\keyword{models}
\keyword{datagen}
