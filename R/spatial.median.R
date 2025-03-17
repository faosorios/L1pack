## ID: spatial.median.R, last updated 2024-05-17, F.Osorio

spatial.median <-
function(x, data, subset, na.action, tol = 1e-6, maxiter = 200)
{ ## generalized spatial median estimator
  ## Rao (1988), Sankhya A 50, 289-313. 
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
    stop("spatial.median applies only to numerical variables")
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
  fit <- .C("spatial_median",
            z = z,
            n = as.integer(n),
            p = as.integer(p),
            median = as.double(center),
            Scatter = as.double(Scatter),
            distances = as.double(distances),
            weights = as.double(rep(1, n)),
            logLik = as.double(logLik),
            tol = as.double(tol),
            maxiter = as.integer(maxiter),
            iterations = as.integer(0))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              x = z,
              dims = dz,
              median = fit$median,
              Scatter = matrix(fit$Scatter, ncol = p),
              logLik = fit$logLik,
              numIter = fit$maxiter,
              innerIter = fit$iterations,
              weights = fit$weights,
              distances = fit$distances,
              speed = speed,
              converged = FALSE)
  names(out$median) <- znames
  dimnames(out$Scatter) <- list(znames, znames)
  if (out$numIter < maxiter)
    out$converged <- TRUE

  class(out) <- "spatial.median"
  out
}

print.spatial.median <-
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
  if (x$converged) {
    out <- character()
    out <- c(out, paste("Converged in", x$numIter, "iterations "))
    out <- c(out, paste("(total inner iterations ", x$innerIter, ")", sep = ""))
    cat(strwrap(paste(out, collapse = "")), sep = "\n")
  } else
    cat("Maximum number of iterations exceeded")
  cat("\nSpatial median:\n ")
  print(format(round(x$median, digits = digits)), quote = F, ...)
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
