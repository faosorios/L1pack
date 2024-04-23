## ID: wilson_hilferty.R, last updated 2024-01-03, F.Osorio

WH.Laplace <- function(x, center, Scatter)
{ # Wilson-Hilferty transformation
  if (missing(center) && missing(Scatter)) {
    if (is.null(x$distances) && is.null(x$dims))
      stop("x is not a valid object.")
    distances <- x$distances
    n <- x$dims[1]
    p <- x$dims[2]
  } else {
    distances <- Mahalanobis(x, center, cov = Scatter, inverted = FALSE)
    n <- nrow(x)
    p <- ncol(x)
    if (p != nrow(Scatter))
      stop("Scatter matrix is not compatible.")
    if (p != length(center))
      stop("center vector is not compatible.")
  }

  z <- .C("Wilson_Hilferty_Laplace",
          distances = as.double(distances),
          n = as.integer(n),
          p = as.integer(p),
          z = double(n))$z
  z
}
