lad <-
function(formula, data, method = c("BR", "EM"), subset, na.action,
  control = NULL, model = TRUE, x = FALSE, y = FALSE, contrasts = NULL)
{
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$method <- mf$control <- mf$model <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)
  ynames <- names(y)
  xnames <- dimnames(x)[[2]]
  dx <- dim(x)
  n  <- dx[1]

  ## set control values
  if (is.null(control))
    control <- l1pack.control()
  method <- match.arg(method)
  choice <- switch(method, "BR" = 0, "EM" = 1)
  ctrl  <- unlist(control)
  nctrl <- names(control)
  ctrl  <- c(ctrl, choice, 0)

  ## initial estimates
  fit <- lsfit(x, y, intercept = FALSE)
  res <- fit$residuals
  cf <- fit$coefficients
  R <- qr.R(fit$qr)

  ## Call fitter
  now <- proc.time()
  fit <- .C("lad",
            y = as.double(y),
            x = as.double(x),
            dims = as.integer(dx),
            coefficients = as.double(cf),
            scale  = as.double(0),
            fitted = double(n),
            resid  = as.double(res),
            weights = as.double(rep(1, n)),
            control = as.double(ctrl),
            sad = as.double(0),
            logLik = as.double(0))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              dims = fit$dims,
              coefficients = fit$coefficients,
              scale = fit$scale,
              minimum = fit$sad,
              fitted.values = fit$fitted,
              residuals = fit$resid,
              R = R,
              numIter = fit$control[4],
              control = fit$control,
              weights = fit$weights,
              logLik = fit$logLik,
              speed = speed,
              converged = FALSE)
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  names(out$control) <- c(nctrl, "method", "numIter")
  names(out$coefficients) <- xnames
  names(out$residuals) <- ynames
  names(out$fitted) <- ynames
  names(out$weights) <- ynames
  dimnames(out$R) <- list(xnames, xnames)
  out$na.action <- attr(mf, "na.action")
  out$contrasts <- attr(x, "contrasts")
  out$xlevels <- .getXlevels(Terms, mf)
  out$terms <- Terms
  if (model)
    out$model <- mf
  if (ret.y)
    out$y <- y
  if (ret.x)
    out$x <- x
  class(out) <- "lad"
  out
}

print.lad <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded")
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$scale), "\n")
  invisible(x)
}

residuals.lad <-
function(object, type = c("working", "response", "quantile"), ...)
{
  type <- match.arg(type)
  r <- object$residuals
  res <- switch(type,
                working =, response = r,
                quantile = {
                  y <- model.response(object$model)
                  mu <- fitted(object)
                  dispersion <- object$scale
                  logp <- plaplace(y, location = mu, scale = dispersion, log.p = TRUE)
                  qnorm(logp, log.p = TRUE)
                }
         )
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
  attr(val, "nobs") <- N ## basic observations must be deleted?
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
  rs <- residuals(x, type = "quantile") ## a kind of standardized residual
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
  if (id.n > 0) { ## label the largest residuals
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
  ##---------- Do the individual plots : ----------
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
  if (show[2L]) { ## Laplace
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

summary.lad <-
function (object, ...)
{
  z <- object
  p <- z$dims[2]
  storage.mode(z$R) <- "double"
  o <- .C("lad_acov",
          R = z$R,
          dims = as.integer(z$dims),
          acov = double(p^2))
  cov.unscaled <- matrix(o$acov, ncol = p)
  se <- sqrt(0.5 * diag(cov.unscaled))
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
  class(ans) <- "summary.lad"
  ans
}

print.summary.lad <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)
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
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$scale))
  cat("\nLog-likelihood:", format(x$logLik), "on", p + 1, "degrees of freedom\n")
  invisible(x)
}
