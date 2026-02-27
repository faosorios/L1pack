\name{residuals.lad}
\alias{residuals.lad}
\title{Residuals for LAD regression}
\description{
  Computes raw or randomized quantile residual for least absolute regression.
}
\usage{
\method{residuals}{lad}(object, type = c("working", "response", "quantile"), ...)
}
\arguments{
  \item{object}{ a fitted model object.}
  \item{type}{ a character string indicating the \code{type} of residual. 
    Available options are \code{"working"}, \code{"response"} for raw residuals, 
    and \code{"quantile"} for quantile residuals, the latter being standardized 
    residuals.}
  \item{\dots}{ additional argument(s) for related methods.}
}
\value{
  Residuals extracted from the \code{lad} object.
}
\references{
  Dunn, P.K., Smyth, G.K. (1996). 
  Randomized quantile residuals. 
  \emph{Journal of Computational and Graphical Statistics} \bold{5}, 236-244. 
  \doi{10.1080/10618600.1996.10474708}
}
\examples{
data(ereturns)
fm <- lad(m.marietta ~ CRSP, data = ereturns)
r <- residuals(fm, type = "quantile")
obs <- c(8,15,34)
plot(r, ylim = c(-2,6), ylab = "quantile residual")
abline(h = 0, lwd = 2, col = "gray75")
abline(h = c(-2,2), lwd = 2, lty = 2, col = "red")
text(obs, r[obs], as.character(obs), pos = 3)
}
\keyword{models}
