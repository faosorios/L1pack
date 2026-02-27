\name{ereturns}
\docType{data}
\alias{ereturns}

\title{Excess returns for Martin Marietta and American Can companies}
\description{
  Data from the Martin Marietta and American Can companies collected over a period 
  of 5 years on a monthly basis.
}
\usage{data(ereturns)}
\format{
  A data frame with 60 observations on the following 4 variables.
  \describe{
    \item{Date}{ the month in which the observations were collected.}
    \item{am.can}{ excess returns from the American Can company.}
    \item{m.marietta}{ excess returns from the Martin Marietta company.}
    \item{CRSP}{ an index for the excess rate returns for the New York stock exchange.}
  }
}
\source{
  Butler, R.J., McDonald, J.B., Nelson, R.D., and White, S.B. (1990).
  Robust and partially adaptive estimation of regression models.
  \emph{The Review of Economics and Statistics} \bold{72}, 321-327. 
  \doi{10.2307/2109722} 
}
\keyword{datasets}