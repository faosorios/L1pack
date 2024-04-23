## ID: distn.R, last updated 2023-12-30, F.Osorio and T.Wolodzko

dlaplace <- function(x, location = 0, scale = 1, log = FALSE)
{ # density of the Laplace distribution
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
{ # distribution function of the Laplace distribution
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
{ # quantile function of the Laplace distribution
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

dmLaplace <-
function(x, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)), log = FALSE)
{ # pdf for the multivariate Laplace distribution
  give.log <- log
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(center))
    stop("'center' must be provided")
  if (isFALSE(center))
    center <- double(p) # center is zeroed
  if (!is.vector(center))
    stop("'center' must be a vector")
  if (length(center) != p)
    stop("'center' has incorrect length")

  if (is.null(Scatter))
    stop("'Scatter' matrix must be provided")
  if (!is.matrix(Scatter))
    stop("'Scatter' must be a matrix")
  if (all(dim(Scatter) != c(p,p)))
    stop("'Scatter' has incorrect dimensions")
  if (!isSymmetric(Scatter))
    Scatter <- asSymmetric(Scatter)

  storage.mode(x) <- "double"
  storage.mode(Scatter) <- "double"

  y <- .C("pdf_mlaplace",
          y = double(n),
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(center),
          Scatter = Scatter)$y

  if (!give.log) y <- exp(y)
  y
}

