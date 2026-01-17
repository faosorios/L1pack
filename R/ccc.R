## ID: ccc.R, last updated 2025-06-26, F.Osorio

l1ccc <- function(x, data, equal.means = FALSE, boots = TRUE, nsamples = 1000, subset, na.action)
{ # estimation of the L1 concordance correlation coefficient
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
    mf$equal.means <- mf$boots <- mf$nsamples <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  }
  else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("'l1ccc' applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  if (p > 2)
    stop("'l1ccc' is not implemented for p > 2")
  
  ## computing L1 and Lin's versions of CCC under Laplace distribution
  o <- Laplace.rho1(z, equal.means)

  ## computing L1 version of CCC under normal distribution
  L1 <- gaussian.rho1(z)

  ## estimation using distribution free U-statistics (King & Chinchilli, 
  ## J. Biopharm. Stat. 11, 83-105, 2001; Stat. Med. 20, 2131-2147, 2001)
  u <- ustat.rho1(z)

  ## computing variances using Bootstrap
  if (boots)
    bvar <- boots.rho1(z, equal.means, nsamples)

  ## creating the output object
  out <- list(call = Call, x = z, dims = dz, rho1 = o$rho1, L1 = list(Laplace = o$rho1,
              Gaussian = L1, U.stat = u$rho1), Lin = o$Lin, ustat = list(rho1 = u$rho1, 
              var.rho1 = u$var.rho1, ustat = u$u, cov = u$v), center = o$center, 
              Scatter = o$Scatter, logLik = o$logLik, weights = o$weights, 
              equal.means = equal.means, boots = boots)
  if (equal.means)
    out$Restricted <- o$Restricted
  if (boots) {
    out$SE <- sqrt(bvar)
    names(bvar) <- NULL
    out$Lin$var.ccc <- bvar[2]
    out$Lin <- out$Lin[c(1,6,2:5)]
    out$Gaussian <- list(rho1 = L1, var.rho1 = bvar[3])
    out$var.rho1 <- bvar[1]
    if (equal.means) {
      out$Restricted$var.ccc <- bvar[4]
      out$Restricted <- out$Restricted[c(1,8,2:7)]
      out <- out[c(1:4,17,5,16,6:11,15,12:14)]
    } else 
      out <- out[c(1:4,16,5,15,6:11,14,12:13)]
  }

  class(out) <- "L1ccc"
  out
}

print.L1ccc <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)

  # L1 coefficients
  if (x$boots) {
    cf <- matrix(c(x$rho1, sqrt(x$var.rho1), 
            x$Gaussian$rho1, sqrt(x$Gaussian$var.rho1), 
            x$ustat$rho1, sqrt(x$ustat$var.rho1)), ncol = 3)
    colnames(cf) <- c("Laplace","Gaussian","U-statistic")
    rownames(cf) <- c("estimate","std.err.")
  } else {
    cf <- c(x$rho1, x$L1$Gaussian, x$L1$U.stat)
    names(cf) <- c("Laplace","Gaussian","U-statistic")
  }
  cat("\nL1 coefficients using:\n ")
  print(format(round(cf, digits = digits)), quote = F, ...)

  # Lin's coefficient
  if (x$boots) {
    cf <- c(x$Lin$ccc, x$Lin$var.ccc, x$Lin$accuracy, x$Lin$precision)
    names(cf) <- c("estimate","variance","accuracy","precision")
  } else {
    cf <- c(x$Lin$ccc, x$Lin$accuracy, x$Lin$precision)
    names(cf) <- c("estimate","accuracy","precision")
  }
  cat("\nLin's coefficients:\n ")
  print(format(round(cf, digits = digits)), quote = F, ...)

  invisible(x)
}

Laplace.rho1 <- function(x, equal.means)
{ ## L1 version of CCC under Laplace distribution
  o <- LaplaceFit(x)
  center <- o$center
  Scatter <- o$Scatter
  cnames <- names(center)
  names(center) <- dimnames(Scatter) <- NULL

  # computing Lin's CCC under Laplace distribution
  p <- 2
  phi  <- Scatter[lower.tri(Scatter, diag = TRUE)]
  diff <- center[1] - center[2]
  ratio <- phi[1] / phi[3]
  rect <- 4 * (p + 1)
  a <- diff / (sqrt(rect) * (phi[1] * phi[3])^.25) # location shift
  b <- sqrt(ratio) # scale shift
  rhoc <- 2 * rect * phi[2] / (rect * (phi[1] + phi[3]) + diff^2)
  accu <- 2 / (b + 1 / b + a^2)
  r    <- phi[2] / sqrt(phi[1] * phi[3]) # correlation coefficient

  # for calculation of the truncated expectation under Laplace distribution
  fnc <- function(x) x * dlaplace(x, scale = 2 * sqrt(2))

  # computing L1 version of CCC under Laplace distribution
  tau <- sqrt(phi[1] + phi[3] - 2 * phi[2])
  ratio <- diff / tau
  ex <- integrate(fnc, lower = -Inf, upper = -ratio)$value
  pr <- plaplace(-ratio, scale = 2 * sqrt(2))
  num <- diff * (1 - 2 * pr) - 2 * tau * ex
  sigma <- sqrt(phi[1] + phi[3])
  ratio <- diff / sigma
  ex <- integrate(fnc, lower = -Inf, upper = -ratio)$value
  pr <- plaplace(-ratio, scale = 2 * sqrt(2))
  den <- diff * (1 - 2 * pr) - 2 * sigma * ex
  rho1 <- 1 - num / den

  # restricted estimation under equality of means
  if (equal.means)
    rf <- restricted.rho1(o)

  ## creating the output object
  names(center) <- cnames
  dimnames(Scatter) <- list(cnames, cnames)
  out <- list(rho1 = rho1, Lin = list(ccc = rhoc, accuracy = accu, 
              precision = r, shift.location = a, shift.scale = b), 
              center = center, Scatter = Scatter, logLik = o$logLik,
              weights = o$weights)
  if (equal.means)
    out$Restricted <- rf
  out
}

gaussian.rho1 <- function(x)
{ ## L1 version of CCC under normal distribution
  o <- cov.weighted(x)
  xbar <- o$mean
  xcov <- o$cov
  names(xbar) <- dimnames(xcov) <- NULL
  phi <- xcov[lower.tri(xcov, diag = TRUE)]
  diff <- xbar[1] - xbar[2]
  tau <- sqrt(phi[1] + phi[3] - 2 * phi[2])
  ratio <- diff / tau
  num <- diff * (1 - 2 * pnorm(-ratio)) + tau * sqrt(2 / pi) * exp(-.5 * ratio^2)
  sigma <- sqrt(phi[1] + phi[3])
  ratio <- diff / sigma
  den <- diff * (1 - 2 * pnorm(-ratio)) + sigma * sqrt(2 / pi) * exp(-.5 * ratio^2)
  rho1 <- 1 - num / den
  rho1
}

restricted.rho1 <- function(x)
{ ## restricted estimation under equality of means and computation of rho1
  call <- x$call
  tol <- 1e-6
  maxiter <- 200
  z <- x$x
  n <- x$dims[1]
  p <- x$dims[2]
  storage.mode(z) <- "double"
  now <- proc.time()
  fit <- .C("fitter_EQUAL",
            x = z,
            n = as.integer(n),
            p = as.integer(p),
            center = as.double(x$center),
            lambda = as.double(mean(x$center)),
            Scatter = as.double(x$Scatter),
            distances = as.double(x$distances),
            weights = as.double(x$weights),
            logLik = as.double(x$logLik),
            tol = as.double(tol),
            maxiter = as.integer(maxiter))
  speed <- proc.time() - now
  # creating the output object for restricted estimation
  Fit <- list(call = call, x = z, dims = x$dims, center = rep(fit$lambda, p),
              lambda = fit$lambda, Scatter = matrix(fit$Scatter, ncol = p),
              weighted.center = fit$center, logLik = fit$logLik, weights = fit$weights, 
              distances = sqrt(fit$distances), numIter = fit$maxiter, speed = speed, 
              converged = FALSE)

  # computing restricted Lin's CCC
  phi  <- Fit$Scatter[lower.tri(Fit$Scatter, diag = TRUE)]
  rho0 <- 2. * phi[2] / (phi[1] + phi[3])
  b0 <- sqrt(phi[1] / phi[3]) # scale shift
  acc0 <- 2. / (b0 + 1 / b0)
  r0 <- phi[2] / sqrt(phi[1] * phi[3])

  # computing asymptotic variance
  rel <- rho0 / r0
  var.rho0 <- 0.25 * (1 - r0^2) * (1 - rho0^2) / (1 - rho0)
  var.rho0 <- var.rho0 * rel^2
  var.rho0 <- var.rho0 / (n - 2)

  # 
  names(Fit$center) <- names(x$center)
  dimnames(Fit$Scatter) <- dimnames(x$Scatter)
  if (Fit$numIter < maxiter)
    Fit$converged <- TRUE
  class(Fit) <- "LaplaceFit"

  ## creating the output object
  out <- list(ccc = rho0, rho1 = 1 - sqrt(1 - rho0), var.rho1 = var.rho0, 
              accuracy = acc0, precision = r0, shifts = list(location = 0, 
              scale = b0), Fitted = Fit)
  out
}

boots.rho1 <- function(x, equal.means, nsamples)
{ ## Bootstrap estimation of standard error
  obs <- 1:nrow(x)
  stat <- matrix(0, nrow = nsamples, ncol = 4)
  cnames <- c("Laplace","Lin","Gaussian","Restricted")

  # initialize progress bar
  cat(" Bootstrap progress:\n")
  pb <- txtProgressBar(min = 0, max = nsamples, style = 3)
  for (i in 1:nsamples) {
    star <- sample(obs, replace = TRUE)
    o <- Laplace.rho1(x[star,], equal.means)
    stat[i,1] <- o$rho1
    stat[i,2] <- o$Lin$ccc
    stat[i,3] <- gaussian.rho1(x[star,])
    if (equal.means)
      stat[i,4] <- o$Restricted$ccc
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  if (!equal.means) {
    cnames <- c("Laplace","Lin","Gaussian")
    stat <- stat[,-4]
  }
  colnames(stat) <- cnames
  out <- apply(stat, 2, var)
  out
}

ustat.rho1 <- function(x)
{ ## estimation and asymptotic variance of rho1 based on U-statistics
  z <- x
  x <- z[,1]
  y <- z[,2]
  n <- nrow(z)

  # computing underlying kernels for the U-statistics
  z <- .Fortran("rho1_ustat",
                x = as.double(x),
                y = as.double(y),
                n = as.integer(n),
                p1 = double(n),
                p2 = double(n))

  # estimating U-statistics and their covariance
  phi <- cbind(z$p1, z$p2)
  z <- cov.weighted(phi)
  u <- z$mean
  v <- (4 / n) * z$cov

  # estimating rho1 using U-statistics
  u1 <- u[1]
  u2 <- u[2]
  h = (n - 1) * (u2 - u1)
  g = u1 + (n - 1) * u2
  rho1 = h / g

  # asymptotic variance of rho1 based on U-statistics
  v11 <- v[1,1]
  v12 <- v[1,2]
  v22 <- v[2,2]
  var.h = (v11 + v22 - 2 * v12) * (n - 1)^2
  var.g = v22 * (n - 1)^2 + v11 + 2 * (n - 1) * v12
  cov.hg = (n - 1) * ((n - 1) * v22 - v11 - (n - 2) * v12)
  var.rho1 = rho1^2 * (var.h / h^2 + var.g / g^2 - 2 * cov.hg / (h * g))

  # output object  
  o <- list(rho1 = rho1, var.rho1 = var.rho1, kernel = phi, u = u, v = v)
  o
}
