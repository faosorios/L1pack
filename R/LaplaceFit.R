## ID: LaplaceFit.R, last updated 2023-12-31, F.Osorio

LaplaceFit <-
function(x, data, subset, na.action, tol = 1e-6, maxiter = 200)
{ ## multivariate-Laplace fitter
  Call <- match.call()
  if (missing(x))
    stop("'x' is not supplied")
  if (inherits(x, "formula")) {
    mt <- terms(x, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "x"] <- "formula"
    mf$tol <- mf$maxiter <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  } else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("LaplaceFit applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  storage.mode(z) <- "double"
 
  ## initial estimates
  o <- cov.weighted(z)
  center <- o$mean
  Scatter <- o$cov
  distances <- Mahalanobis(z, center, Scatter)
  R <- chol(Scatter)
  logLik <- -0.5 * n * p * log(2 * pi) - n * sum(log(diag(R))) - 0.5 * sum(distances)

  ## Call fitter
  now <- proc.time()
  fit <- .C("Laplace_fitter",
            z = z,
            n = as.integer(n),
            p = as.integer(p),
            center = as.double(center),
            Scatter = as.double(Scatter),
            distances = as.double(distances),
            weights = as.double(rep(1, n)),
            logLik = as.double(logLik),
            tol = as.double(tol),
            maxiter = as.integer(maxiter))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              x = z,
              dims = dz,
              center = fit$center,
              Scatter = matrix(fit$Scatter, ncol = p),
              logLik = fit$logLik,
              numIter = fit$maxiter,
              weights = fit$weights,
              distances = sqrt(fit$distances),
              speed = speed,
              converged = FALSE)
  names(out$center) <- znames
  dimnames(out$Scatter) <- list(znames, znames)
  if (out$numIter < maxiter)
    out$converged <- TRUE

  class(out) <- "LaplaceFit"
  out
}

print.LaplaceFit <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <- function(z, digits = digits, ...) {
    # print upper triangular part of covariance matrix
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded")
  cat("\nCenter:\n ")
  print(format(round(x$center, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  if (x$dims[2] <= 5)
    print.symmetric(x$Scatter, digits = digits)
  else {
    print.symmetric(x$Scatter[1:5,1:5], digits = digits)
    cat("...")
  }
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs, "\n")
  invisible(x)
}
