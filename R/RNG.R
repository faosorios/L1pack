## ID: RNG.R, last updated 2023-05-17, F.Osorio and T.Wolodzko

rlaplace <- function(n, location = 0, scale = 1)
{ # univariate Laplace random generation
  if (scale <= 0.0)
    stop("'scale' must be non-negative.")

  nloc <- length(location)
  nscale <- length(scale)
  # calling C code
  x <- .C("r_laplace",
          n = as.integer(n),
          x = double(n),
          location = as.double(location),
          nloc = as.integer(nloc),
          scale = as.double(scale),
          nscale = as.integer(nscale))$x
  x
}

rmLaplace <-
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
{ # multivariate Laplace random number generation
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # calling C code
  y <- .C("RNG_mlaplace",
          y = as.double(y),
          dims = as.integer(dy),
          center = as.double(center),
          Scatter = as.double(Scatter))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}
