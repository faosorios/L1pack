## ID: l1fit.R, last updated 2022-10-18, F.Osorio

l1fit <-
function(x, y, intercept = TRUE, tolerance = 1e-07, print.it = TRUE)
{ ## wrapper to Barrodale and Roberts procedure
  ## based on S-Plus 'l1fit' function.
  xn <- if (is.matrix(x)) dimnames(x)[[2]] else "X"
  x <- as.matrix(x)
  dx <- dim(x)
  vars <- dx[2]
  if (length(xn) == 0)
    xn <- paste("X", 1:vars, sep = "")
  if (intercept) {
    vars <- vars + 1
    xn <- c("Intercept", xn)
    x <- cbind(1, x)
  }
  nobs <- dx[1]
  if (nobs != length(y))
    stop("x and y have different number of observations")
  if (nobs <= vars)
    stop("More variables than observations")
  storage.mode(y) <- "double"
  z <- matrix(0, nrow = nobs + 2, ncol = vars + 2)
  z[1:nobs, 1:vars] <- x
  storage.mode(z) <- "double"
  fit <- .Fortran("l1",
                  n = as.integer(nobs),
                  p = as.integer(vars),
                  n2 = as.integer(nobs + 2),
                  p2 = as.integer(vars + 2),
                  z = z,
                  y = y,
                  tol = as.double(tolerance),
                  cf = double(vars),
                  resid = double(nobs),
                  work = integer(nobs))
  minimum <- fit$z[nobs + 1, vars + 1]
  rank <- fit$z[nobs + 1, vars + 2]
  info <- fit$z[nobs + 2, vars + 1]
  iter <- fit$z[nobs + 2, vars + 2]
  fitted <- drop(x %*% fit$cf)
  msg <- character(0)
  if (info == 0) {
    msg <- "Non-unique solution possible."
    if (print.it)
      warning(msg)
  }
  else if(info != 1)
    stop("Premature termination.")
  if (rank != vars) {
    msg <- c(msg, paste("Matrix not of full rank, apparent rank", rank))
    if (print.it)
      warning(msg[length(msg)])
  }

  # creating the output object
  out <- list(coefficients = fit$cf,
              minimum = minimum,
              fitted.values = fitted,
              residuals = fit$resid,
              rank = rank,
              iterations = iter,
              info = info)
  names(out$coefficients) <- xn
  if (length(msg))
    out$message <- list(msg)
  out
}
