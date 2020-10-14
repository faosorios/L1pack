## ID: distn.R, last updated 2020-10-09, F.Osorio and T.Wolodzko

dlaplace <- function(x, location = 0, scale = 1, log = FALSE)
{
  if (scale <= 0.0)
    stop("'scale' must be non-negative.")

  n <- length(x)
  nloc <- length(location)
  nscale <- length(scale)

  y <- .C("d_laplace",
          n = as.integer(n),
          y = double(n),
          x = as.double(x),
          location = as.double(location),
          nloc = as.integer(nloc),
          scale = as.double(scale),
          nscale = as.integer(nscale),
          give.log = as.integer(log))$y
  y
}

plaplace <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
  if (scale <= 0.0)
    stop("'scale' must be non-negative.")

  n <- length(q)
  nloc <- length(location)
  nscale <- length(scale)

  y <- .C("p_laplace",
          n = as.integer(n),
          y = double(n),
          q = as.double(q),
          location = as.double(location),
          nloc = as.integer(nloc),
          scale = as.double(scale),
          nscale = as.integer(nscale),
          lower.tail = as.integer(lower.tail),
          log.p = as.integer(log.p))$y
  y
}

qlaplace <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
  if (scale <= 0.0)
    stop("'scale' must be non-negative.")

  n <- length(p)
  nloc <- length(location)
  nscale <- length(scale)

  y <- .C("q_laplace",
          n = as.integer(n),
          y = double(n),
          p = as.double(p),
          location = as.double(location),
          nloc = as.integer(nloc),
          scale = as.double(scale),
          nscale = as.integer(nscale),
          lower.tail = as.integer(lower.tail),
          log.p = as.integer(log.p))$y
  y
}

rlaplace <- function(n, location = 0, scale = 1)
{
  if (scale <= 0.0)
    stop("'scale' must be non-negative.")

  nloc <- length(location)
  nscale <- length(scale)

  x <- .C("r_laplace",
          n = as.integer(n),
          x = double(n),
          location = as.double(location),
          nloc = as.integer(nloc),
          scale = as.double(scale),
          nscale = as.integer(nscale))$x
  x
}
