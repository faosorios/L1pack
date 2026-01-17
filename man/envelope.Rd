\name{envelope.Laplace}
\alias{envelope.Laplace}
\title{QQ-plot with simulated envelopes}
\description{
Constructs a normal QQ-plot using a Wilson-Hilferty transformation for the estimated Mahalanobis distances obtained from the \code{\link{LaplaceFit}} procedure.
}
\usage{
envelope.Laplace(object, reps = 50, conf = 0.95, plot.it = TRUE)
}
\arguments{
\item{object}{ an object of class \code{\link{LaplaceFit}} representing the fitted model.}
\item{reps}{ number of simulated point patterns to be generated when computing the envelopes. The default number is 50.}
\item{conf}{ the confidence level of the envelopes required. The default is to find 95\% confidence envelopes.}
\item{plot.it}{ if \code{TRUE} it will draw the corresponding plot, if \code{FALSE} it will only return the computed values. }
}
\value{
  A list with the following components :
  
  \item{transformed}{ a vector with the \code{z}-scores obtained from the Wilson-Hilferty transformation.}
  \item{envelope}{ a matrix with two columns corresponding to the values of the lower and upper pointwise 
  confidence envelope.}
}
\references{
  Atkinson, A.C. (1985).
  \emph{Plots, Transformations and Regression}.
  Oxford University Press, Oxford.

  Osorio, F., Galea, M., Henriquez, C., Arellano-Valle, R. (2023). 
  Addressing non-normality in multivariate analysis using the t-distribution. 
  \emph{AStA Advances in Statistical Analysis} \bold{107} 785-813.

  Vallejos, R., Osorio, F., Ferrer, C. (2025+).
  A new coefficient to measure agreement between continuous variables.
  Working paper.
}
\keyword{dplot}
