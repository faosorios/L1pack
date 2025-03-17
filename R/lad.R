## ID: lad.R, last updated 2024-01-16, F.Osorio and T.Wolodzko

lad <-
function(formula, data, subset, na.action, method = "BR", tol = 1e-7, maxiter = 200, 
  x = FALSE, y = FALSE, contrasts = NULL)
{ ## least absolute deviations fit
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$method <- mf$tol <- mf$maxiter <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)

  # call fitter
  z <- lad.fit(x, y, method, tol, maxiter)

  # output
  z$call <- Call
  z$model <- mf
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(Terms, mf)
  z$terms <- Terms
  if (ret.y)
    z$y <- y
  if (ret.x)
    z$x <- x
  class(z) <- "lad"
  z
}

lad.fit <-
function(x, y, method = "BR", tol = 1e-7, maxiter = 200)
{ ## dispatcher among various fitting functions
  if (!is.numeric(x))
    stop("model matrix must be numeric.")
  if (!is.numeric(y))
    stop("response must be numeric.")
  if (!length(x))
    method <- "null"
  contr <- attr(x, 'contrast')
  fit <- switch(method,
                BR  = lad.fit.BR(x, y, tol),
                EM  = lad.fit.EM(x, y, tol, maxiter),
                stop(paste("unimplemented method:", method)))
  fit$contrast <- contr
  fit$method <- method
  fit
}

lad.fit.BR <-
function(x, y, tol = 1e-7)
{
  if (is.matrix(y))
    stop("'lad.fit.BR' does not support multiple responses.")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
	n <- dx[1]
  ynames <- names(y)
	xnames <- dimnames(x)[[2]]
	if (n != ny)
		stop("Number of observations in x and y not equal.")
	p <- dx[2]
  if (n <= p)
    stop("More variables than observations.")

  # call fitter
  converged <- TRUE
  now <- proc.time()
  z <- .C("lad_BR",
          y = y,
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          coef = double(p),
          scale = as.double(0),
          sad = as.double(0),
          fitted = double(n),
          residuals = double(n),
          logLik = as.double(0),
          tol = as.double(tol),
          rank = as.integer(0),
          iter = as.integer(0),
          info = as.integer(0))
  speed <- proc.time() - now
  msg <- character(0)
  if (z$info == 0)
    msg <- "Non-unique solution possible."
  else if(z$info != 1) {
    msg <- "Premature termination."
    converged <- FALSE
  }
  if (z$rank != p)
    msg <- c(msg, paste("Matrix not of full rank, apparent rank", z$rank))

  # post-processing
  basic <- (1:n)[z$residuals == 0.0]
  R <- chol(crossprod(x))

  # creating the output object
  out <- list(coefficients = z$coef, scale = z$scale, residuals = z$residuals,
              fitted.values = z$fitted, SAD = z$sad, logLik = z$logLik, basic = basic,
              dims = dx, R = R, iterations = z$iter, speed = speed, converged = converged)
  if (length(msg))
    out$message <- msg
  names(out$coefficients) <- xnames
  names(out$residuals) <- ynames
  names(out$fitted.values) <- ynames
  dimnames(out$R) <- list(xnames, xnames)
  class(out) <- "lad"
  out
}

lad.fit.EM <-
function(x, y, tol = 1e-7, maxiter = 200)
{
  if (is.matrix(y))
    stop("'lad.fit.EM' does not support multiple responses.")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
	n <- dx[1]
  ynames <- names(y)
	xnames <- dimnames(x)[[2]]
	if (n != ny)
		stop("Number of observations in x and y not equal.")
	p <- dx[2]
  if (n <= p)
    stop("More variables than observations.")

  # set control values
  control <- c(maxiter, tol, 0)

  # initial estimates
  fit <- lsfit(x, y, intercept = FALSE)
  res <- fit$residuals
  cf <- fit$coefficients
  R <- qr.R(fit$qr)

  # call fitter
  converged <- FALSE
  now <- proc.time()
  z <- .C("lad_EM",
          y = y,
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          coefficients = as.double(cf),
          scale = as.double(0),
          fitted = double(n),
          residuals = as.double(res),
          weights = as.double(rep(1, n)),
          sad = as.double(0),
          logLik = as.double(0),
          tol = as.double(tol),
          iter = as.integer(maxiter))
  speed <- proc.time() - now

  # post-processing
  if (z$iter < maxiter)
    converged <- TRUE
  eps <- .Machine$double.eps^.35
  dev <- z$residuals
  test <- abs(dev) < eps * z$sad
  weights <- ifelse(test, 1.0, eps / abs(dev))
  basic <- (1:n)[test]
  msg <- character(0)
  if (length(basic) == 0) {
    basic <- NULL
    msg <- "method was unable to determine basic observations"
  }

  # creating the output object
  out <- list(coefficients = z$coef, scale = z$scale, residuals = z$residuals,
              fitted.values = z$fitted, SAD = z$sad, logLik = z$logLik,
              weights = weights, basic = basic, dims = dx, R = R, iterations = z$iter,
              speed = speed, converged = converged)
  if (length(msg))
    out$message <- msg
  names(out$coefficients) <- xnames
  names(out$residuals) <- ynames
  names(out$fitted.values) <- ynames
  names(out$weights) <- ynames
  dimnames(out$R) <- list(xnames, xnames)
  class(out) <- "lad"
  out
}

print.lad <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call)
  if (x$converged)
    cat("Converged in", x$iterations, "iterations\n")
  else {
    if (x$method == "EM")
      cat("Maximum number of iterations exceeded\n")
  }
  cat("\nCoefficients:\n")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$scale), "\n")
  if (!is.null(x$message))
    cat(paste("\nNOTE:", x$message, "\n"))
  invisible(x)
}

residuals.lad <-
function(object, type = c("working", "response", "quantile"), ...)
{
  type <- match.arg(type)
  r <- object$residuals
  res <- switch(type,
                working = r, response = r,
                quantile = {
                  y <- model.response(object$model)
                  mu <- fitted(object)
                  dispersion <- object$scale
                  logp <- plaplace(y, location = mu, scale = dispersion, log.p = TRUE)
                  qnorm(logp, log.p = TRUE)
                })
  res
}

fitted.lad <- function(object, ...)
  object$fitted.values

## sqrt(2) 'factor' has not been included
deviance.lad <- function(object, ...)
  sum(abs(residuals(object)), na.rm = TRUE)

## log-likelihood for lad objects
logLik.lad <- function(object, ...)
{
  val <- object$logLik
  N <- object$dims[1]
  p <- object$dims[2]
  attr(val, "nall") <- N
  attr(val, "nobs") <- N # basic observations must be deleted?
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}

predict.lad <- function(object, newdata, na.action = na.pass, ...)
{
  if (missing(newdata) || is.null(newdata))
    return(fitted(object))
  tt <- terms(object)
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)
  X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  drop(X %*% coef(object))
}

plot.lad <-
function(x, which = 1:2,
         caption = c("Residuals vs Fitted", "Laplace Q-Q plot"),
         panel = points,
         sub.caption = deparse(x$call), main = "",
         ask = prod(par("mfcol")) < length(which) && dev.interactive(),
         ...,
         id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75)
{
  if (!inherits(x, "lad"))
    stop("Use only with 'lad' objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 2))
    stop("`which' must be in 1:2")
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  rs <- residuals(x, type = "quantile") # a kind of standardized residual
  yhat <- fitted(x)
  rs.lab <- "Standardized residuals"
  yh.lab <- "Fitted values"
  n <- length(rs)
  if (is.null(id.n))
    id.n <- 0
  else {
    id.n <- as.integer(id.n)
    if (id.n < 0 || id.n > n)
      stop("`id.n' must be in {1,..,",n,"}")
  }
  if (id.n > 0) { # label the largest residuals
    if (is.null(labels.id))
      labels.id <- paste(1:n)
    iid <- 1:id.n
    show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
    text.id <- function(x, y, ind, adj.x = FALSE)
      text(x - if(adj.x) strwidth(" ")*cex.id else 0, y, labels.id[ind],
           cex = cex.id, xpd = TRUE, adj = if(adj.x) 1)
  }
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  #---------- Do the individual plots : ----------
  if (show[1L]) {
    ylim <- range(rs, na.rm = TRUE)
    if (id.n > 0)
      ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)
    plot(yhat, rs, xlab = yh.lab, ylab = rs.lab, main = main, ylim = ylim, type = "n", ...)
    panel(yhat, rs, ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[1], 3, 0.25)
    if (id.n > 0) {
      y.id <- rs[show.rs]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      text.id(yhat[show.rs], y.id, show.rs, adj.x = TRUE)
    }
    abline(h = 0, lty = 3, col = "gray")
  }
  if (show[2L]) { # Laplace
    ylim <- range(rs, na.rm=TRUE)
    ylim[2] <- ylim[2] + diff(ylim) * 0.075
    qq <- qqnorm(rs, main = main, ylab = rs.lab, ylim = ylim, ...)
    abline(c(0,1), lty = 3, col = "gray")
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[2], 3, 0.25)
    if (id.n > 0)
      text.id(qq$x[show.rs], qq$y[show.rs], show.rs, adj.x = TRUE)
  }
  if (!one.fig && par("oma")[3] >= 1)
  mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}

simulate.lad <- function(object, nsim = 1, seed = NULL, ...)
{
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)  # initialize the RNG if necessary
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ftd <- fitted(object)
  nm <- names(ftd)
  n <- length(ftd)
  ntot <- n * nsim
  ans <- as.data.frame(ftd +
                       matrix(rlaplace(ntot, scale = object$scale), nrow = n))
  names(ans) <- paste("sim", seq_len(nsim), sep = "_")
  if (!is.null(nm))
    row.names(ans) <- nm
  attr(ans, "seed") <- RNGstate
  ans
}

confint.lad <- function(object, parm, level = 0.95, ...)
{
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) parm <- seq(along = pnames)
  else if(is.character(parm)) parm <- match(parm, pnames, nomatch = 0)
  a <- (1 - level) / 2
  a <- c(a, 1-a)
  pct <- paste(round(100 * a, 1), "%")
  ci <- array(NA, dim = c(length(parm), 2),
              dimnames = list(pnames[parm], pct))
  se <- sqrt(diag(vcov(object)))[parm]
  qz <- qnorm(a)
  ci[] <- cf[parm] + se %o% qz
  ci
}

vcov.lad <- function(object, ...)
{
  scale <- object$scale
  2.0 * scale^2 * chol2inv(object$R)
}

summary.lad <-
function (object, ...)
{
  z <- object
  p <- z$dims[2]
  cov.unscaled <- chol2inv(z$R)
  se <- z$scale * sqrt(2.0 * diag(cov.unscaled))
  est <- z$coefficients
  zval <- est / se
  ans <- z[c("call", "terms")]
  ans$dims <- z$dims
  ans$logLik <- z$logLik
  ans$scale <- z$scale
  ans$residuals <- z$residuals
  ans$coefficients <- cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients),
        c("Estimate", "Std.Error", "Z value", "p-value"))
  ans$cov.unscaled <- cov.unscaled
  class(ans) <- "summary.lad"
  ans
}

print.summary.lad <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call)
  resid <- x$residuals
  nobs <- x$dims[1]
  p <- x$dims[2]
  rdf <- nobs - p
  if (rdf > 5) {
    cat("\nResiduals:\n")
		rq <- quantile(resid)
		names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
		print(rq, digits = digits, ...)
	}
	else if(rdf > 0) {
	 cat("\nResiduals:\n")
	 print(resid, digits = digits, ...)
  }
  cat("\nCoefficients:\n")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$scale))
  cat("\nLog-likelihood:", format(x$logLik), "on", p + 1, "degrees of freedom\n")
  invisible(x)
}
